// gravity engine by mushy
#include <allegro.h>
#include <math.h>
#define W 1024
#define H 768
#define MODE GFX_AUTODETECT_FULLSCREEN
//#define INITIAL_VELOCITY_ON 5
//#define SUN_ON TRUE
//#define BOUNCE ON
//#define MOMENTUM_EXCHANGE
//#define MINIGRAVITY
//#define POLAR_FORCES
//#define GRAVITY
//#define STRONG_FORCE
#define SHM_GRAVITY
#define NUM 200
#define S 0.00000000001
double scale = S/NUM*NUM;
// ball radius, x, y, and y velocity
struct particle
{
	double x;
	double y;
	double xv;
	double yv;
	double m;
	int rd;
	int r, g, b;
	double col;
	int particle_type;
}p[NUM];
double attract(struct particle *a, struct particle *b)
{
	// scale * mass / dist
	double mag=scale*b->m;
#ifdef POLAR_FORCES
	if(a->particle_type == b->particle_type)
	{
		mag *= -1;
	}
#endif

	///(((a->x - b->x) * (a->x - b->x)) + ((a->y - b->y) * (a->y - b->y)));
	double distsq = ((a->x - b->x) * (a->x - b->x)) + ((a->y - b->y) * (a->y - b->y));
	double dist = sqrt(distsq);
	if (a->rd + b->rd < dist){
	//mag*=(b->m*b->m*b->m*b->m);
		double magx = mag * abs(a->x - b->x)
#ifdef SHM_GRAVITY
								* dist
#endif
								/(dist
#ifdef GRAVITY
				*distsq
#endif
#ifdef MINIGRAVITY
				*dist
#endif
				);


		double magy = mag * abs(a->y - b->y)
#ifdef SHM_GRAVITY
								* dist
#endif
								/(dist
#ifdef GRAVITY
				*distsq
#endif
#ifdef MINIGRAVITY
				*dist
#endif
				);

#ifdef STRONG_FORCE
		magx += mag * abs(a->x - b->x) * 10000000 / (distsq * distsq) ;
		magy += mag * abs(a->y - b->y) * 10000000 / (distsq * distsq) ;
#endif

		a->xv -= a->x > b->x ? magx : -magx;
		a->yv -= a->y > b->y ? magy : -magy;
	} else {
#ifdef BOUNCE
#	ifdef MOMENTUM_EXCHANGE
		double saved_momentum=(b->xv * b->m)/a->m;
		b->xv=(a->xv * a->m)/b->m;
		a->xv=saved_momentum;
		saved_momentum=(b->yv * b->m)/a->m;
		b->yv=(a->yv * a->m)/b->m;
		a->yv=saved_momentum;
#	else
		a->xv = -a->xv;
		a->yv = -a->yv;
#	endif
#endif
	}
	a->x += a->xv;
	a->y += a->yv;
#ifdef POTENTIAL_ON
	return dist * scale * b->m * a->m/NUM;
#endif
}
// double buffer
BITMAP *buffer;
int main(void)
{
	allegro_init();
	install_timer();
	install_mouse();
	install_keyboard();
	srand(time(NULL));
	set_gfx_mode(MODE, W, H, 0, 0);
	buffer = create_bitmap(W, H);
	int i, j;
	double total_mass=0;
	double centre_of_mass_x, centre_of_mass_y, last_cog_x, last_cog_y;
	double white=makecol(255,255,255);
	double total_kinetic, total_potential;
	double highscore, pot_highscore,start_potential;
	for(i = 0; i < NUM; i++)
	{
		p[i].x = rand() % W;
		p[i].y = rand() % H;
#ifdef INITIAL_VELOCITY_ON
		p[i].xv = ((rand() % 100) - 100) / 100;
		p[i].yv = ((rand() % 100) - 100) / 100;
#endif
		p[i].r = rand() % 255;
		p[i].g = rand() % 255;
		p[i].b = rand() % 255;
		p[i].rd = (rand() % 15) + 1;
		p[i].m = p[i].rd*p[i].rd*p[i].rd;//*p[i].rd*p[i].rd*p[i].rd;
		p[i].col=makecol(p[i].r, p[i].g, p[i].b);
		p[i].particle_type = rand() % 2;
		total_mass+= p[i].m;
	}
#ifdef SUN_ON
	p[0].x=W/2;
	p[0].y=H/2;
	p[0].rd = 20;
	p[0].m = 50000;
	p[0].col=white;
	p[0].xv=0;
	p[0].yv=0;
#endif

	highscore = 0;
	pot_highscore= 0;
	last_cog_x = W/2;
	last_cog_y = H/2;
	while(!key[KEY_ESC])
	{
		clear_bitmap(buffer);
		total_mass = 0;
		centre_of_mass_x=0;
		centre_of_mass_y=0;
		total_kinetic=0;
		total_potential=0;
		for(i = 0; i < NUM; i++)
		{
			for(j = 0; j < NUM; j++)
			{
				if(i != j)
				{
					total_potential += attract(&p[i],&p[j]);
				}
			}
			if(i==0)
			{
				total_mass = p[0].m;
				centre_of_mass_x=p[0].x;
				centre_of_mass_y=p[0].y;
			}
			else
			{
				total_mass += p[i].m;
				centre_of_mass_x+=(p[i].m/total_mass)*(p[i].x-centre_of_mass_x);
				centre_of_mass_y+=(p[i].m/total_mass)*(p[i].y-centre_of_mass_y);
			}
			total_kinetic += (p[i].xv * p[i].xv + p[i].yv * p[i].yv) * p[i].m/2;
			circlefill(buffer, p[i].x, p[i].y, p[i].rd, p[i].col);
		}

		highscore = MAX(highscore, total_kinetic);
#ifdef POTENTIAL_ON
		pot_highscore = MAX(total_potential, pot_highscore);
		if (start_potential < 1){
			start_potential = total_potential;
		}
#endif
//		last_cog_x = centre_of_mass_x;
//		last_cog_y = centre_of_mass_y;

		textprintf(buffer, font, 0, 0, 15, "KINETIC ENERGY: %f", total_kinetic);
		textprintf(buffer, font, 0, 10, 15, "HIGHSCORE: %f", highscore);
#ifdef POTENTIAL_ON
		textprintf(buffer, font, 0, 20, 15, "POTENTIAL: %f", total_potential);
		textprintf(buffer, font, 0, 30, 15, "POT HIGHSCORE: %f", pot_highscore);
		textprintf(buffer, font, 0, 40, 15, "TOTAL ENERGY: %f", total_kinetic + total_potential);
		textprintf(buffer, font, 0, 50, 15, "start_potential: %f", start_potential);
#endif

		circle(buffer, centre_of_mass_x,centre_of_mass_y, 20, white );
		blit(buffer, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
		//rest(10);
	}
	allegro_exit();
}
END_OF_MAIN()
