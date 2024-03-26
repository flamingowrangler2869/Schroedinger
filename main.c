#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

#define WIDTH	512		// width of the window
#define HEIGHT	300		// height of the window
#define PSISCL	10		// position y-axis scaling factor
#define PHISCL	0.1		// same but momentum
#define RES		256		// resolution
#define COEFFCT	64		// amount of coeffs

#define M	1	// mass of the particle
#define L	1	// width of square well
#define H	1	// planck's constant

const double pi = acos(-1);		// pi
const double hbar = H/(2*pi);	// reduced planck's constant

double	t = 0.0;	// elapsed time
int32_t mouse_x;	// mouse x position
int32_t mouse_y;	// mouse y position
bool	mouse_down;	// if mouse button 1 is held
bool	wf_view = true;	// real+imag or prob+phase
size_t	nth_ef = 0;		// view the nth eigenfunction
bool	paused = false;	// simulation paused
double	speed = 0.05;	// simulation speed

complex double psi[RES];		// wave function
complex double phi[RES];		// momentum (fourier transform of wavefunction)
complex double coeffs[COEFFCT]; // superpositional coefficients
double states[COEFFCT];			// energy states E_n
SDL_Renderer *renderer;

#define SQR(x)		((x)*(x))				// square of x
#define FLIP(a) 	((a)=!(a))				// opposite of a
#define LERP(s,a,b)	((1-(s))*(a)+((s)*(b))) // linear interpolation
#define MAX(a,b)	((a)>(b)?(a):(b))		// maximum of two values
#define MIN(a,b)	((a)<(b)?(a):(b))		// minimum of two values

double eigenfunction(size_t n, int x) {
	return sqrt(2.0/L) * sin(n*pi*x/(L*RES));
}

void fft(complex double in[], complex double out[], int step, int n) {
	/* in:   time domain		unaffected by transform
	 * out:  frequency domain	output
	 * step: slicing mechanism	always 1 in foreign calls
	 * n:    length of array	must be power of 2
	 */
	if (n == 1) {
		out[0] = in[0];
		return;
	}
	fft(in,		 out, 	  step*2, n/2);
	fft(in+step, out+n/2, step*2, n/2);
	for (int k = 0; k < n/2; k++) {
		float t = (float)k/n;
		complex double w = cexp(-2*I*pi*t) * out[k+n/2];
		complex double e = out[k];
		out[k]     = e + w;
		out[k+n/2] = e - w;
	}
}

double minval(complex double *array, size_t len) {
	double current_min = 0;
	for (size_t i = 0; i < len; i++) {
		double n;
		if ((n=creal(array[i])) < current_min) current_min = n;
		if ((n=cimag(array[i])) < current_min) current_min = n;
	}
	return current_min;
}

double maxval(complex double *array, size_t len) {
	double current_max = 0;
	for (size_t i = 0; i < len; i++) {
		double n;
		if ((n=creal(array[i])) > current_max) current_max = n;
		if ((n=cimag(array[i])) > current_max) current_max = n;
	}
	return current_max;
}

void clear_func(complex double array[], size_t len) {
	for (size_t i = 0; i < len; i++) array[i] = 0.0;
}

void swap_halves(complex double array[], size_t len) {
	complex double tmp[len/2];
	for (size_t i = 0; i < len/2; i++) {
		tmp[i] = array[i];
	}
	for (size_t i = 0; i < len; i++) {
		if (i < len/2)
			array[i] = array[i+len/2];
		else
			array[i] = tmp[i-len/2];
	}
}

void set_color(uint32_t color) {
	/* color = 0x**RRGGBB */
	uint8_t r, g, b;
	r = (color>>16)&0xFF;
	g = (color>> 8)&0xFF;
	b = (color>> 0)&0xFF;
	SDL_SetRenderDrawColor(renderer, r, g, b, 255);
}

double map(double x, double i_min, double i_max, double o_min, double o_max) {
	return o_min + (o_max-o_min)/(i_max-i_min)*(x-i_min);
}

double map0(double x, double i_max, double o_max) {
	return map(x, 0, i_max, 0, o_max);
}

uint32_t hue_to_rgb(size_t hue) {
	/* hue: 0..359
	 * return: 0x00RRGGBB
	 */
	hue %= 360;
	float hue_pr = hue/60.0;
	float x = 255*(1-fabsf(fmodf(hue_pr, 2.0)-1));
	struct Color { uint8_t r, g, b; } color1 = {0, 0, 0};
	if		(hue_pr >= 0 && hue_pr < 1) color1 = (struct Color){255, x, 0};
	else if (hue_pr >= 1 && hue_pr < 2)	color1 = (struct Color){x, 255, 0};
	else if (hue_pr >= 2 && hue_pr < 3)	color1 = (struct Color){0, 255, x};
	else if (hue_pr >= 3 && hue_pr < 4)	color1 = (struct Color){0, x, 255};
	else if (hue_pr >= 4 && hue_pr < 5)	color1 = (struct Color){x, 0, 255};
	else if (hue_pr >= 5 && hue_pr < 6) color1 = (struct Color){255, 0, x};
	return (uint32_t)(((color1.r)<<16) | ((color1.g)<<8) | (color1.b));
}

double probability(complex double wavefunc[], size_t x) {
	return SQR(cabs(wavefunc[x]));
}

void draw_probability(complex double *func, float scale, int axis, int ceiling) {
	axis = HEIGHT/2-axis;
	for (int x = 0; x < WIDTH; x++) {
		size_t n0 = (size_t)map0(x,  WIDTH,RES);
		size_t n1 = (size_t)map0(x+1,WIDTH,RES);
		double phase = carg(func[n0])+pi; // 0 rad = 1+0i
		double prob0 = probability(psi, n0);
		double prob1 = probability(psi, n1);
		double interpolated = LERP((map0(x,WIDTH,RES)-n0), prob0, prob1);
		double height = interpolated;
		set_color(hue_to_rgb(map0(phase,pi,180)));
		SDL_RenderDrawLine(renderer,
			x, axis,
			x, axis-(int)(scale*height)
		);
	}
}

void draw_complex(complex double *func, float scale, int axis, int ceiling, int color_r, int color_i) {
	axis = HEIGHT/2-axis;
	set_color(color_i);
	for (int x = 0; x < WIDTH; x++) {
		SDL_RenderDrawLine(renderer,
			    x*(float)(WIDTH/(RES-1)), axis+(int)(scale*cimag(func[x  ])),
			(x+1)*(float)(WIDTH/(RES-1)), axis+(int)(scale*cimag(func[x+1]))
		);
	}
	set_color(color_r);
	for (int x = 0; x < RES; x++) {
		SDL_RenderDrawLine(renderer,
			    x*(float)(WIDTH/(RES-1)), axis+(int)(scale*creal(func[x  ])),
			(x+1)*(float)(WIDTH/(RES-1)), axis+(int)(scale*creal(func[x+1]))
		);
	}
}

double expectation(complex double wavefunc[], size_t len) {
	double sum = 0.0;
	for (size_t x = 0; x < len; x++)
		sum += x * probability(wavefunc,x) / len;
	return sum;
}

#ifdef WINDOWS
int WinMain() {
#else
int main() {
#endif
	SDL_Init(SDL_INIT_VIDEO);
	SDL_Window   *window   = SDL_CreateWindow("Particle in a Box", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL);
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_Event event;
	const Uint8 *keyboard_state = SDL_GetKeyboardState(NULL);

reset:
	// initial conditions
	t = 0.0;
	double amplitude = HEIGHT/25.0;
	double spread = (mouse_y == 0) ? 0.01 : exp(mouse_y*0.02);
	for (int x = 0; x < RES; x++) {
		psi[x] = -amplitude*exp(-SQR(((x-mouse_x/2)/spread)));
	}

	if (!mouse_x && !mouse_y) {
		// program start
		coeffs[1] = 5;
		coeffs[2] = 5;
	} else {
		for (size_t n = 0; n < COEFFCT; n++) {
			coeffs[n] = 0.0;
			for (size_t x = 0; x < RES; x++)
				coeffs[n] += psi[x] * conj(eigenfunction(n,x)) / RES;
//			// Eₙ = (n²π²ħ²) / (2mL²)
			states[n] = SQR(n*pi*hbar)/(2*M*SQR(L));
		}
	}

	Uint32 then = SDL_GetTicks();
	for (size_t frame = 0 ;; !paused ? frame++ : 0) {
		SDL_GetMouseState(&mouse_x, &mouse_y);
		SDL_PollEvent(&event);
		if (keyboard_state[SDL_SCANCODE_ESCAPE] || event.type == SDL_QUIT) goto exit;
		switch (event.type) {
			case SDL_MOUSEBUTTONDOWN:
			case SDL_MOUSEBUTTONUP: mouse_down = event.button.state; break;
			case SDL_KEYDOWN:
				switch (event.key.keysym.sym) {
					case SDLK_UP:	speed *= 1.03; break;
					case SDLK_DOWN:	speed /= 1.03; break;
					case SDLK_SPACE: FLIP(paused); break;
					case SDLK_r: goto reset;
					case SDLK_w: FLIP(wf_view); break;
				} break;
		}
		set_color(0);
		SDL_RenderClear(renderer);

		for (size_t x = 0; x < RES; x++) {
			psi[x] = 0.0;
			for (size_t n = 0; n < COEFFCT; n++) {
				psi[x] += coeffs[n] * eigenfunction(n, x) * cexp(-I*states[n]*t*speed/hbar);
			}
		}

		fft(psi, phi, 1, RES);
		swap_halves(phi, RES);

		if (wf_view)
			draw_complex(psi, PSISCL,  0, 50/*ceiling*/, 0xFFFFFF, 0x0000FF);
		else
			draw_probability(psi, 1, -HEIGHT/2/*-HEIGHT/4*/, 0/*ceiling*/ );
		int expect = (int)expectation(psi, RES);
		set_color(0x00FFFFFF);
//		SDL_RenderDrawLine(renderer, expect, 0, expect, HEIGHT);
		SDL_RenderPresent(renderer);

#if 1
		printf("frame\t%zu %s"		"\e[K\n"
			   "time\t%.3lfs"		"\e[K\n"
			   "speed\t%.2lfx"		"\e[K\n"
			   "expect\t%d"			"\e[K\n"
			   "mouse\t(%d, %d) %s"	"\e[K"
			   "\r\e[4A",
			frame, paused ? "PAUSED" : "",
			t*speed, speed, expect,
			mouse_x, mouse_y, mouse_down ? "down" : "up"
		);
#endif

		if (mouse_down) goto reset;

		// timing
		Uint32 now = SDL_GetTicks();
		double dt = now - then;
		if (!paused) t += dt/1000.0;
		then = now;
	}

exit:
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
	printf("\n\n\n\n\n");
	return 0;
}
