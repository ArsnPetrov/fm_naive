#include <SDL2/SDL.h>
#include <rtl-sdr.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <fftw3.h>

#define SAMPLE_RATE 2048000
#define BUFFER_LENGTH 16384 * 4
#define SAMPLE_RATE_AUDIO 48000

#define MHz(x) x * 1000000
#define kHz(x) x * 1000

int filter = 0;

void IQ_to_int8(unsigned char *buf, int8_t *res, uint32_t len)
{
	for (int i = 0; i < len; ++i)
	{
		res[i] = (int)buf[i] - 128;
	}
}

int8_t get(int8_t *array, int32_t i)
{
	if (i < 0)
	{
		return 0;
	}
	else 
	{
		return array[i];
	}
}

float get(float *array, int32_t i, float *memory)
{
	if (i < 0)
	{
		return memory[-i - 1];
	}
	else 
	{
		return array[i];
	}
}

const double lpf_80kHz_1_a[] = { 1 };

void lowpass80(int8_t *iq, int8_t *res, uint32_t len)
{
	int order;
	double coeff_a[9] = { 0 };
	double coeff_b[9] = { 0 };
	long double a_denom;
	//double a_denom = 307644801.280;
	// double a_denom = 681.663;
	if (filter == 0)
	{
		a_denom = 9.108;
		order = 1;
		coeff_a[0] = 1;
		coeff_b[0] = 0.7804128239;
	}
	else if (filter == 1) {
		a_denom = 681.663;
		coeff_a[0] = 3;
		coeff_a[1] = 3;
		coeff_a[2] = 1;
		coeff_b[0] = 1711.182 / a_denom;
		coeff_b[1] = -1454.237 / a_denom;
		coeff_b[2] = 416.718 / a_denom;
		order = 3;
	}
	else if (filter == 2) {
		a_denom = 307644801.280;
		coeff_a[0] = 9;
		coeff_a[1] = 36;
		coeff_a[2] = 84;
		coeff_a[3] = 126;
		coeff_a[4] = 126;
		coeff_a[5] = 84;
		coeff_a[6] = 36;
		coeff_a[7] = 9;
		coeff_a[8] = 1;
		coeff_b[0] = 2334063651.004 / a_denom;
		coeff_b[1] = -7899938578.743 / a_denom;
		coeff_b[2] = 15652560522.784 / a_denom;
		coeff_b[3] = -20003674508.411 / a_denom;
		coeff_b[4] = 17096789287.654 / a_denom;
		coeff_b[5] = -9770813786.735 / a_denom;
		coeff_b[6] =  3599991706.286 / a_denom;
		coeff_b[7] = -775839999.813 / a_denom;
		coeff_b[8] = 74505995.254 / a_denom;
		order = 8;
	}

	
	// double coeff_a[] = {
	// 	9, 
	// 	36,
	// 	84,
	// 	126,
	// 	126,
	// 	84,
	// 	36,
	// 	9, 
	// 	1
	// };
	// double coeff_b[] = {
	// 	7.5868782482, 
	// 	-25.6787650754, 
	// 	50.8786771551, 
	// 	-65.0219812758,
	// 	55.5731454473,
	// 	-31.7600484262,
	// 	11.7017797515,
	// 	-2.5218693655,
	// 	0.2421818764
	// };
	double _I, _Q;
	for (int i = 0; i < len / 2; ++i)
	{
		_I = (double)get(iq, 2*i);
		_Q = (double)get(iq, 2*i);
		for (int j = 1; j <= order; j++)
		{
			_I += coeff_a[j - 1] * get(iq, 2*(i - j));
		}
		_I /= a_denom;
		for (int j = 1; j <= order; j++)
		{
			_I += coeff_b[j - 1] * get(res, 2*(i - j));
		}

		for (int j = 1; j <= order; j++)
		{
			_Q += coeff_a[j - 1] * get(iq, 2*(i - j) + 1);

		}
		_Q /= a_denom;
		for (int j = 1; j <= order; j++)
		{
			_Q += coeff_b[j - 1] * get(res, 2*(i - j) + 1);
		}

		res[2 * i] = (int8_t)_I;
		res[2 * i + 1] = (int8_t)_Q;
	}
}

void phase_diff(int8_t *iq, float *res, uint32_t len)
{
	float a, b, c, d, re, im;

	for (int i = 1; i < len / 2; ++i)
	{
		a = iq[2 * i];
		b = -iq[2 * i + 1]; // сопряжённое

		c = iq[2 * (i + 1)];
		d = iq[2 * (i + 1) + 1];

		re = a * c - b * d;
		im = a * d + b * c;

		res[i] = atan2f(im, re);
	}

	res[0] = res[1];
}

void decimate(float *a, uint32_t len, float *res, int factor)
{
	for (int i = 0; i < len / factor; ++i)
	{
		res[i] = a[i * factor];
	}
}

void upsample(float *a, uint32_t len, float *res, int factor)
{
	memset(res, 0, len * factor);

	for (int i = 0; i < len; i++)
	{
		res[i * factor] = a[i];
	}
}

void lowpass_256(float *data, float *res, uint32_t len, float *memory)
{
	for (int i = 0; i < len; ++i)
	{
		//res[i] = get(data, i);
		// res[i] = (get(data, i, memory)
		// + 3 * get(data, i-1, memory)
		// + 3 * get(data, i-2, memory)
		// + 1 * get(data, i-2, memory)
		// + 2122805.632 * get(res, i-1, memory)
		// + (-2091204.823) * get(res, i-1, memory)
		// + 686767.715 * get(res, i-2, memory)) / 718376.525;

		res[i] = (get(data, i, memory)
		+ 1 * get(data, i-1, memory)
		+ 87.892 * get(res, i-2, memory)) / 89.892;
	}

	// memory[2] = res[len - 3];
	// memory[1] = res[len - 2];
	// memory[0] = res[len - 1];
}

void divide_a_by_n(float *a, uint32_t len, float n)
{
	for (int i = 0; i < len; ++i)
	{
		a[i] /= n;
	}
}

struct demod_ctx {
	int received_buffer_length;
	int audio_buffer_length;

	bool *quit;
	bool updated;
	bool audio_updated;

	int8_t received_buffer   [BUFFER_LENGTH];
	int8_t lowpassed         [BUFFER_LENGTH];
	float audio_long         [BUFFER_LENGTH / 2];
	float audio_interpolated [BUFFER_LENGTH / 2 * 3];
	float audio_lowpassed    [BUFFER_LENGTH / 2 * 3];
	float audio_decimated    [BUFFER_LENGTH / 2 * 3 / 128];
	float audio              [BUFFER_LENGTH / 2 * 3 / 128];

	fftw_complex audio_long_complex[BUFFER_LENGTH / 2];
	fftw_complex audio_long_fft[BUFFER_LENGTH / 2];
	
	fftw_complex audio_complex[BUFFER_LENGTH / 2 * 3 / 128];
	fftw_complex audio_fft[BUFFER_LENGTH / 2 * 3 / 128];

	double last_rms;
	double last_mean;

	rtlsdr_dev_t *device;
	pthread_t *device_thread;

	uint32_t frequency;
};

void intHandler(int a)
{
	exit(-1);
}

void device_callback(unsigned char *buf, uint32_t len, void *ctx)
{
	//printf("sdr got a new buffer\n");

	struct demod_ctx *_ctx = (struct demod_ctx*)ctx;

	_ctx->updated = true;

	IQ_to_int8(buf, ((demod_ctx*)ctx)->received_buffer, len);
}

static void *device_thread_fn(void *arg)
{
	struct demod_ctx *_ctx = (struct demod_ctx*)arg;
	rtlsdr_dev_t *_dev = _ctx->device;

	rtlsdr_read_async(_dev, device_callback, arg, 0, BUFFER_LENGTH);

	return 0;
}

static void *freq_reader_fn(void *arg)
{
	struct demod_ctx *_ctx = (struct demod_ctx*)arg;
	float f = 0;
	while(1)
	{
		scanf("%f", &f);
		_ctx->frequency = (int)kHz(f);
		rtlsdr_cancel_async(_ctx->device);
		rtlsdr_set_center_freq(_ctx->device, _ctx->frequency);
		rtlsdr_reset_buffer(_ctx->device);
		rtlsdr_set_tuner_gain_mode(_ctx->device, 0);
		printf("f = %d\n", _ctx->frequency);
		pthread_create(_ctx->device_thread, NULL, device_thread_fn, (void*)(_ctx));
	}
}

static void *demod_thread_fn(void *arg)
{
	struct demod_ctx *_ctx = (struct demod_ctx*)arg;

	while (1) {
		if (_ctx->updated) {
			_ctx->updated = false;

			double rms = 0;
			double mean = 0;
			double mad = 0;

			lowpass80(_ctx->received_buffer, _ctx->lowpassed, BUFFER_LENGTH);
			phase_diff(_ctx->lowpassed, _ctx->audio_long, BUFFER_LENGTH);
			upsample(_ctx->audio_long, BUFFER_LENGTH / 2, _ctx->audio_interpolated, 3);
			float memory[3] = { 0 };
			lowpass_256(_ctx->audio_interpolated, _ctx->audio_lowpassed, BUFFER_LENGTH / 2 * 3, memory);
			decimate(_ctx->audio_lowpassed, BUFFER_LENGTH / 2 * 3, _ctx->audio, 128);

			for (int i = 0; i < BUFFER_LENGTH / 2 * 3 / 128; ++i)
			{
				rms += powf(_ctx->audio[i], 2) / (BUFFER_LENGTH / 2 * 3 / 128);
				mean += _ctx->audio[i] / (BUFFER_LENGTH / 2 * 3 / 128);
				mad += fabs(_ctx->audio[i]) / (BUFFER_LENGTH / 2 * 3 / 128);
			}

			rms = sqrt(rms);

			_ctx->last_rms = rms;
			_ctx->last_mean = mean;

			divide_a_by_n(_ctx->audio, BUFFER_LENGTH / 2 * 3 / 128, 2 * mad);

			_ctx->audio_updated = true;
		}
	}

	return 0;
}

void audio_device_callback(void *param, Uint8 *stream, int len)
{
	float *s = (float*)stream;
	float *audio = (float*)param;
	
	memcpy(stream, audio, len);
}

int main()
{
	signal(SIGINT, intHandler);

	int r;
	rtlsdr_dev_t *device;
	pthread_t device_thread, demodulator_thread, frequency_reader_thread;
	bool quit = false;
	struct demod_ctx ctx;
	SDL_AudioSpec spec;
	char *audio_device_name;
	int audio_device_id = 1;
	SDL_Event e;
	SDL_Window *window;
	SDL_Renderer *renderer;
	fftw_plan audio_fft_plan;

	SDL_Init(SDL_INIT_EVERYTHING);

	SDL_CreateWindowAndRenderer(2000, 600, SDL_WINDOW_SHOWN | SDL_WINDOW_ALLOW_HIGHDPI, &window, &renderer);

	spec.freq = SAMPLE_RATE_AUDIO;
    spec.format = AUDIO_F32SYS;
    spec.channels = 1;
    spec.samples = BUFFER_LENGTH / 2 * 3 / 128;
    spec.callback = audio_device_callback;
    spec.userdata = ctx.audio;

    audio_device_name = (char*)SDL_GetAudioDeviceName(audio_device_id, 0);

    auto recorder = SDL_OpenAudioDevice(audio_device_name, 0, &spec, &spec, SDL_AUDIO_ALLOW_FORMAT_CHANGE);

    if (recorder == 0) {
        printf("[fm-receiver] Could not connect to the audio device: %s\n", SDL_GetError());
        return 1;
    }

    // 

	ctx.updated = false;
	ctx.audio_updated = false;
	ctx.received_buffer_length = BUFFER_LENGTH;
	ctx.audio_buffer_length = 48000;
	ctx.frequency = MHz(103.4);
	ctx.device_thread = &device_thread;

	audio_fft_plan = fftw_plan_dft_1d(BUFFER_LENGTH / 2 * 3 / 128, ctx.audio_complex, ctx.audio_fft, FFTW_FORWARD, FFTW_ESTIMATE);

	r = rtlsdr_open(&device, 0);

	if (r)
	{
		printf("[fm-receiver] Could not open RTL-SDR. Error code: %d\n", r);
		exit(-1);
	}

	rtlsdr_reset_buffer(device);
	rtlsdr_set_sample_rate(device, SAMPLE_RATE);
	rtlsdr_set_center_freq(device, ctx.frequency);
	rtlsdr_set_tuner_gain_mode(device, 0);
	printf("[fm-receiver] Central frequency is set to %d kHz\n", rtlsdr_get_center_freq(device) / 1000);

	ctx.device = device;
	
	pthread_create(&device_thread, NULL, device_thread_fn, (void*)(&ctx));
	pthread_create(&demodulator_thread, NULL, demod_thread_fn, (void*)(&ctx));
	pthread_create(&frequency_reader_thread, NULL, freq_reader_fn, (void*)(&ctx));
	usleep(100000);
	SDL_PauseAudioDevice(recorder, 0);

	double average_fft[BUFFER_LENGTH / 2 * 3 / 128] = { 0 };

	while(!quit)
	{
		//printf("main thread is alive\n");
		while(SDL_PollEvent(&e))
		{
			if (e.type == SDL_QUIT)
			{
				printf("[fm-receiver] Quitting the program\n");
				rtlsdr_cancel_async(device);
				quit = true;
			}

			if (e.type == SDL_KEYDOWN)
			{
				switch (e.key.keysym.scancode)
				{
				case SDL_SCANCODE_F:
					printf("[fm-receiver] Current frequency is %d Hz.\n", rtlsdr_get_center_freq(device));
					break;

				case SDL_SCANCODE_R:
					printf("[fm-receiver] Last RMS is %f\n", ctx.last_rms);
					printf("[fm-receiver] Last mean is %f\n", ctx.last_mean);
					break;

				case SDL_SCANCODE_0:
					rtlsdr_cancel_async(device);
					rtlsdr_set_center_freq(device, MHz(90.3));
					rtlsdr_reset_buffer(device);
					pthread_create(&device_thread, NULL, device_thread_fn, (void*)(&ctx));
					usleep(100000);
					break;

				case SDL_SCANCODE_1:
					filter++;
					filter %= 3;
					printf("filter set to %d\n", filter);
					break;

				
				default:

					break;
				}
			}
		}

		if (ctx.audio_updated)
		{
			SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
			SDL_RenderClear(renderer);

			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			for (size_t i = 1; i < spec.samples; i++)
			{
				SDL_RenderDrawLine(renderer, i - 1, 200 - 50 * ctx.audio[i - 1], i, 200 -  50 * ctx.audio[i]);
			}
			SDL_RenderDrawLine(renderer, 0, 200, 2000, 200);

			for (size_t i = 1; i < 400; i++)
			{
				SDL_RenderDrawLine(renderer, i - 1, 600 - 1 * ctx.received_buffer[2 * (i - 1)], i, 600 - 1 * ctx.received_buffer[2 * i]);
			}

			for (size_t i = 1; i < 400; i++)
			{
				SDL_RenderDrawLine(renderer, i - 1, 1000 - 1 * ctx.lowpassed[2 * (i - 1)], i, 1000 - 1 * ctx.lowpassed[2 * i]);
			}

			// do fft

			for (size_t i = 0; i < BUFFER_LENGTH / 2 * 3 / 128; i++)
			{
				ctx.audio_complex[i][0] = ctx.audio[i];
				ctx.audio_complex[i][1] = 0;
			}
			
			fftw_execute(audio_fft_plan);

			for (size_t i = 0; i < BUFFER_LENGTH / 2 * 3 / 128; i++)
			{
				// average_fft[i-1] *= 0.95;
				// average_fft[i-1] += fabs(ctx.audio_fft[(i - 1)][0]) * 0.05;
				average_fft[i] *= 0.95;
				average_fft[i] += fabs(ctx.audio_fft[i][0]) * 0.05;
				
				SDL_RenderDrawLine(renderer, 400 + 2*i, 600, 400 + 2*i, 600 - 8 * average_fft[i]);
				// SDL_RenderDrawLine(renderer, 400 + i - 1, 600 - 0.05* average_fft[i-1], 400 + i, 600 - 0.05 * average_fft[i]);

			}

			SDL_RenderDrawLine(renderer, 400, 600, 1400, 600);

			//
			
			SDL_RenderPresent(renderer);

			ctx.audio_updated = false;
		}
	}

	fftw_destroy_plan(audio_fft_plan);
	SDL_PauseAudioDevice(audio_device_id, 1);
	rtlsdr_close(device);

	return 0;
}
