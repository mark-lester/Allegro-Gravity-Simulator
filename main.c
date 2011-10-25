// gravity engine by mushy
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#define ALLEGRO_ON
//#define SDL_ON
//#define NO_GRAPHICS
#ifdef ALLEGRO_ON
#include <allegro.h>
#endif

#ifdef SDL_ON
#include <SDL/SDL.h>
#include <SDL/SDL_draw.h>
#endif

#include <math.h>

#define W 1680
#define H 1050
#define D 32
#define SCREEN_WIDTH W
#define SCREEN_HEIGHT H
#define SCREEN_DEPTH D
#define MODE GFX_AUTODETECT_FULLSCREEN
#define RADIUS_OF_CENTRE_OF_MASS 20
double screen_diag=1;
int world_type=0;
int max_its=-1;
int sun_flag = FALSE;
int bounce=FALSE;
int merge=FALSE;
int initial_velocity_flag=FALSE;
double initial_velocity_range=100;
double extra_factor=0;

#ifdef ALLEGRO_ON
// double buffer
BITMAP *buffer;
#endif

#ifdef SDL_ON
SDL_Surface *screen;

Uint8       *pixel;
#endif


//#define INITIAL_VELOCITY_ON 5
//#define SUN_ON TRUE
//#define BOUNCE ON
//#define MOMENTUM_EXCHANGE
//#define MINIGRAVITY
//#define POLAR_FORCES
//#define GRAVITY
//#define STRONG_FORCE
//#define SHM_GRAVITY


#define POTENTIAL_ON
#define MAX_ARRAY_SIZE 10000
// SHM S 0.000000001,  scale = S/log(NUM)
#define NUM 200
#define S 200
double average_energy_per_object=S;
double scale = 1;
double array_size=200;

// energy = integral of half distance * number of objects
// ball radius, x, y, and y velocity
struct particle
{
	double x;
	double y;
	double xv;
	double yv;
	double dxv;
	double dyv;
	double m;
	int rd;
	int r, g, b;
	double col;
	int particle_type;
}p[MAX_ARRAY_SIZE];

#define FORCE_GRAVITY 0
#define FORCE_MINIGRAVITY 1
#define FORCE_LOG_GRAVITY 2
#define FORCE_CONSTANT 3
#define FORCE_SHMGRAVITY 4
#define FORCE_SUPER_GRAVITY 5
#define FORCE_SUPER3_GRAVITY 6

double merge_it(struct particle *a, struct particle *b){
	if (a->m < b->m){
		struct particle *t;
		t=a;
		a=b;
		b=a;
	}
	a->xv =(b->xv*b->m + a->xv*a->xv)/(a->m + b->m);
	a->yv =(b->yv*b->m + a->yv*a->yv)/(a->m + b->m);

	a->m +=b->m;
	a->rd=sqrt(a->rd*a->rd + b->rd*b->rd);
	b->m=b->rd=0;
	return a->m;
}
double bounce_it(struct particle *a, struct particle *b, double min_dis){
   double rel_vol = sqrt(
		                  ((a->xv - b->xv)*(a->xv - b->xv)) +
		                  ((a->yv - b->yv)*(a->yv - b->yv))
		                  );
#define BSCALE 1
   b->dxv += 2*(a->xv - b->xv)*abs(a->xv - b->xv)/(((b->m/a->m) + 1) * rel_vol * BSCALE);
   b->dyv += 2*(a->yv - b->yv)*abs(a->yv - b->yv)/(((b->m/a->m) + 1) * rel_vol* BSCALE);
   a->dxv += 2*(b->xv - a->xv)*abs(a->xv - b->xv)/(((a->m/b->m) + 1) * rel_vol* BSCALE);
   a->dyv += 2*(b->yv - a->yv)*abs(a->yv - b->yv)/(((a->m/b->m) + 1) * rel_vol* BSCALE);
   double mx = abs(a->xv - b->xv)*min_dis/ (rel_vol*2);
   double my = abs(a->yv - b->yv)*min_dis/ (rel_vol*2);
   mx++;
   my++;
   if (a->x > b->x){
	   a->x += mx;
	   b->x -= mx;
   } else {
	   a->x -= mx;
	   b->x += mx;
   }

   if (a->y > b->y){
	   a->y += my;
	   b->y -= my;
   } else {
	   a->y -= my;
	   b->y += my;
   }

   return rel_vol;
}

double attract(struct particle *a, struct particle *b,int force_type){
	double force,pot;
	double dist = sqrt(((a->x - b->x) * (a->x - b->x)) + ((a->y - b->y) * (a->y - b->y)));
	// no change in velocity if objects overlap to avoid 1/0, and take potential therefore as zero
	double combined_radii = (a->rd + b->rd);
    if ((a->m == 0) || (b->m == 0)){
    	return 0;
    }

	if (combined_radii > dist){
		if (bounce){
			bounce_it(a,b,combined_radii - dist +1);
		}
		if (merge){
			merge_it(a,b);
		}
		return 0;
	}
	//the force  and the potential energy are both proportional to the product of the two masses
	force = a->m * b->m;
	force *= scale;
	pot = force;
	//force *= scale;
    //pot = force;

	// we now need to work out the actual force based on the rules, and integrate over the distance between them to get the potential energy
	switch (force_type){
		case FORCE_CONSTANT:
			// no change for the force "constant", it's constant!
			// force *= 1;

			// integral of a constant is x, i.e. the distance
			pot *= dist - combined_radii;
			break;

		case FORCE_GRAVITY:
			// inverse square
			force /=  dist*dist;

			// integral of 1/(x*x) is -1/x
			pot *= (dist - combined_radii)/(dist*combined_radii);
			break;

		case FORCE_MINIGRAVITY:
			// our mini gravity is just inverse, not inverse square
			force /= dist;

			// integral of 1/x is log(x)
			pot *= log(dist) - log(combined_radii);
			break;

		case FORCE_SHMGRAVITY:
			// SHM gravity is directly prop to distance
			force *= dist /*- combined_radii */;

			// integral of x is x squared over 2
			pot *= (dist*dist  - combined_radii*combined_radii) /2;
			break;

		case FORCE_SUPER_GRAVITY:
			// Super gravity is directly prop to square of distance
			force *= dist * dist;

			// integral of x squared is x cubed over 3
			pot *= (dist*dist*dist - combined_radii*combined_radii*combined_radii) /3;
			break;

		case FORCE_LOG_GRAVITY:
			// log gravity is directly prop to log of distance
			force *= log(dist);

			// integral of log x is x * log(x) - x, or so I am told
			pot *= dist * log(dist) - dist - (combined_radii* log(combined_radii) - combined_radii);
			break;

		case FORCE_SUPER3_GRAVITY:
			force *= dist*dist*dist;

			// integral of x cubed is x to the 4 over 4
			pot *= (dist*dist*dist*dist - combined_radii*combined_radii*combined_radii*combined_radii)/4;
			break;
	}

	// scale the incremental velocity vector components, (diff in x and y comps divided by the distance, e.g. sin and cos)
	// then divide by the mass of the object being moved (newtons 2nd law)
	// do for both a and b
	a->dxv += (force * (b->x - a->x))/(dist * a->m);
	a->dyv += (force * (b->y - a->y))/(dist * a->m);
	b->dxv += (force * (a->x - b->x))/(dist * b->m);
	b->dyv += (force * (a->y - b->y))/(dist * b->m);

	return pot;  // this is the total potential for both objects combined
}

double scale_up(double val){
		if (val < 10){
			val += 1;
		} else {
			val *= 1.1;
		}
		return val;
	}

double scale_down(double val){
		if (val < 10 && val > 2){
			val -= 1;
		} else {
			val /= 1.1;
		}
		return val;
}



int do_world(int force_type, int world_size)
{
	int i, j;
	double total_mass=0;
	double centre_of_mass_x, centre_of_mass_y;
#ifndef NO_GRAPHICS
	double white=makecol(255,255,255);
	double red         = makecol( 255, 0,   0   );
	double green       = makecol( 0,   255, 0   );
//	double blue        = makecol( 0,   0,   255 );
	double yellow      = makecol( 255, 255, 0   );
#endif

	double total_kinetic, total_potential;
	double highscore, pot_highscore,start_potential,start_kinetic,total_energy;
	start_potential=start_kinetic=total_energy=0;
	char *pstring;
	int loop_counter=0;

	int k_finish,p_finish;
	for(i = 0; i < world_size; i++)
	{
		p[i].x = (rand() % (W/2)) + (W/4);
		p[i].y = (rand() % (H/2)) + (H/4);
		p[i].xv = p[i].dxv=  p[i].yv = p[i].dyv=0;

#define RAND_COLOUR (rand() % 200) + 55;
		p[i].r = RAND_COLOUR;
		p[i].g = RAND_COLOUR;
		p[i].b = RAND_COLOUR;
		p[i].rd = (rand() % 15) + 1;
		//p[i].rd = (int)(rand() / 15);
		p[i].m = p[i].rd*p[i].rd*p[i].rd;//*p[i].rd*p[i].rd*p[i].rd;
#ifndef NO_GRAPHICS
		p[i].col=makecol(p[i].r, p[i].g, p[i].b);
#endif
		p[i].particle_type = rand() % 2;
	    if (sun_flag && (i == world_size -1) ){ // make it the last one so it renders on top
	    	p[i].x=W/2;
	    	p[i].y=H/2;
	    	p[i].rd = 25;
	    	p[i].m = 50000;
	    	p[i].col=white;
	    	p[i].xv=0;
	    	p[i].yv=0;
	    }

		total_mass+= p[i].m;
	}


	scale = 1;
	// do a dummy run to add up the potentsail energy at start

	for(i = 0; i < world_size; i++) {
			for(j = i +1 ; j < world_size; j++) {
				  start_potential += attract(&p[i],&p[j],force_type);
			}
			p[i].xv = p[i].dxv=  p[i].yv = p[i].dyv=0; //zap these back to zero
    }

	// this should result in a constant average energy per particle for all models
	scale =  average_energy_per_object*array_size/start_potential;
	start_potential=average_energy_per_object*array_size;

	if (initial_velocity_flag){
		initial_velocity_range=3; // start_potential / total_mass;
		//initial_velocity_range /= (double)10;
		// do a dummy run to add up the potential energy at start
		for(i = 0; i < world_size-2; i+=2) { // dont do the last if it's a sun
			int k,j;
				if (p[i].m < p[i+1].m){
					k=i;j=i+1;
				} else {
					k=i+1;j=i;
				}
			    p[k].xv = (rand() % (int)initial_velocity_range) * (rand() %2 ? 1 : -1) ;
				p[k].yv = (rand() % (int)initial_velocity_range) * (rand() %2 ? 1 : -1);
				p[j].xv = -p[k].xv * p[k].m / p[j].m; // make the next one balance up the momentum
				p[j].yv = -p[i].yv * p[k].m / p[j].m;
				start_kinetic += (p[i].xv * p[i].xv) + (p[i].yv * p[i].yv) * p[i].m ;
				start_kinetic += (p[i+1].xv * p[i+1].xv) + (p[i+1].yv * p[i+1].yv) * p[i+1].m ;
		}
	}


	switch (force_type){
			case FORCE_CONSTANT:
				pstring="constant force = independent of distance ";
				break;

			case FORCE_GRAVITY:
				pstring="normal = inverse square of the distance ";
				break;

			case FORCE_MINIGRAVITY:
				pstring="mini = inverse to the distance  ";
				break;

			case FORCE_SHMGRAVITY:
				pstring="simple harmonic = distance ";
				break;

			case FORCE_SUPER_GRAVITY:
				pstring="super = square of the distance ";
				break;

			case FORCE_LOG_GRAVITY:
				pstring="log = log of the distance";
				break;

			case FORCE_SUPER3_GRAVITY:
				pstring="super 3 = cube of the distance between objects";
				break;
		}


	highscore = 0;
	pot_highscore= 0;
#ifdef ALLEGRO_ON
	while(!key[KEY_ESC] && !key[KEY_M] && !key[KEY_B] && !key[KEY_S] && !key[KEY_V] && !key[KEY_SPACE] && !key[KEY_RIGHT] && !key[KEY_LEFT] && !key[KEY_UP]  && !key[KEY_DOWN]  && !key[KEY_PGUP]  && !key[KEY_PGDN]) {
#endif
#ifdef SDL_ON
	SDL_Event keyevent;    //The SDL event that we will poll to get events.

	while (SDL_PollEvent(&keyevent)){ //Poll our SDL key event for any keystrokes.
	   switch(keyevent.type){
		 case SDL_KEYDOWN:
		    return 0;
	   }
#endif
#ifdef NO_GRAPHICS
	 while (1) {
#endif
		loop_counter++;
		if (max_its > 0){
			if (loop_counter > max_its){
				world_type=rand() %6;
				return 1;
			}
		}
#ifdef ALLEGRO_ON
		while (key[KEY_TAB]){
			rest(10);
		}
		clear_bitmap(buffer);
#endif

		total_mass = 0;
		centre_of_mass_x=0;
		centre_of_mass_y=0;
		total_kinetic=0;
		total_potential=0;
		for(i = 0; i < world_size && p[i].m; i++) {
			// only do this once, so for j > i. attract will work out the new incremental velocity vectors for both objects and return the total potential energy
			for(j = i +1 ; j < world_size; j++) {
				  total_potential += attract(&p[i],&p[j],force_type);
			}
		}

		// we now add on the new velocity components and move the dots about and add up the new kinetic,
		// the potential and kinetic will be one loop out, but it works good enuf
		for(i = 0; i < world_size && p[i].m; i++){
				p[i].xv += p[i].dxv;
				p[i].yv += p[i].dyv;

				p[i].x += p[i].xv;
				p[i].y += p[i].yv;
#ifdef ALLEGRO_ON
				circlefill(buffer, p[i].x, p[i].y, p[i].rd, p[i].col);
#endif
#ifdef SDL_ON
				Draw_FillCircle(screen,
				                     (Sint16)p[i].x, (Sint16) p[i].y, (Uint16) p[i].r,
				                     (Uint32) p[i].col);
#endif

				p[i].dxv = p[i].dyv = 0; // set the incremental velocity vec to 0 for next time

				// work out centre of mass
				if(i==0) {  // first one obviously IS the CoM
					total_mass = p[0].m;
					centre_of_mass_x=p[0].x;
					centre_of_mass_y=p[0].y;
				}
				else {  // then just work out how far you need to shift in the direction of the new body
						// based on what fraction of the total mass the new body is
					total_mass += p[i].m;
					centre_of_mass_x+=(p[i].m/total_mass)*(p[i].x-centre_of_mass_x);
					centre_of_mass_y+=(p[i].m/total_mass)*(p[i].y-centre_of_mass_y);
				}

				// add up kinetic energy
				total_kinetic += ((p[i].xv * p[i].xv) + (p[i].yv * p[i].yv)) * p[i].m/2;

		}

		highscore = MAX(highscore, total_kinetic);
#ifdef POTENTIAL_ON
		pot_highscore = MAX(total_potential, pot_highscore);
		total_energy = total_potential + total_kinetic;

#endif
#ifdef ALLEGRO_ON
		textprintf(buffer, font, 0, 0, 15, "%s: o %d :l %d :v %0.0f: sp %0.0f: sk %0.0f: c %0.0f: k %0.0f: p %0.0f", pstring,world_size,loop_counter, initial_velocity_range,start_potential, start_kinetic, total_energy, total_kinetic, total_potential);

		// total_kinetic_stripe = (total_kinetic / start_potential) * 3/4s of the screen width, same for potential,

		//  start potential energy should equal total energy at all times of course as we started off with no kinetic)
        double start_energy=start_potential+start_kinetic;
		int ypos = 16;
		rectfill(buffer, 0, ypos, W*3/4, ypos+2, green);
		ypos+=3;
		k_finish = total_kinetic*(double)W*3/(4*start_energy);
		p_finish = total_potential*(double)W*3/(4*start_energy);
		p_finish += k_finish + 1;  // draw the potential on the end of the kinetic
		// kinetic energy
		rectfill(buffer, 0, ypos, k_finish, ypos+4, red);
		// potential energy

		rectfill(buffer, k_finish + 1,ypos, p_finish,ypos+4, yellow);

		// draw centre of mass, should stay the same
		circle(buffer, centre_of_mass_x,centre_of_mass_y, RADIUS_OF_CENTRE_OF_MASS, white );
		blit(buffer, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
#endif
#ifdef SDL_ON
		Draw_Circle(screen,
		                     (Sint16)centre_of_mass_x, (Sint16) centre_of_mass_y, (Uint16) RADIUS_OF_CENTRE_OF_MASS,
		                     (Uint32) white);



		/* update the screen (aka double buffering) */
		SDL_Flip(screen);
#endif
		//rest(10);
	}
	 if (loop_counter==0){
		 return 1;
	 }
	 if (key[KEY_B]){
		   bounce = !bounce;
	 }
	 if (key[KEY_M]){
		   merge = !merge;
	 }
	 if (key[KEY_RIGHT]){
	 		 world_type++;
	 }
	 if (key[KEY_LEFT]){
	 		 world_type--;
	 }
	 if (key[KEY_UP]){
		 average_energy_per_object=scale_up(average_energy_per_object);
	 	 }
	 if (key[KEY_DOWN]){
		 average_energy_per_object=scale_down(average_energy_per_object);
	 }

	 if (key[KEY_PGUP]){
	 		 array_size = scale_up(array_size);
	 }
	 if (key[KEY_PGDN]){
		 array_size = scale_down(array_size);
	 }
	 if (key[KEY_S]){
		 sun_flag = !sun_flag;
     }
	 if (key[KEY_V]){
		 initial_velocity_flag = !initial_velocity_flag;
     }

	return !key[KEY_ESC];
}




int main(int argc, char *argv[]){
	int c;

	while ((c = getopt (argc, argv, "mbsn:g:e:i:v:")) != -1) {
	         switch (c)
	           {
	           case 's':
	        	 sun_flag=TRUE;
	        	 break;

	           case 'b':
	        	 bounce=TRUE;
	        	 break;

	           case 'm':
	        	 merge=TRUE;
	        	 break;

	           case 'v':
	        	 initial_velocity_flag=TRUE;
	        	 //initial_velocity_range = atoi(optarg);
	        	 break;

	           case 'n': // number of objects
		         array_size=atoi(optarg);
		         break;

	           case 'g': // gravity type
		         world_type=atoi(optarg);
		         break;

	           case 'e': // average energy
		         average_energy_per_object=atoi(optarg);
	             break;

	           case 'i': // average energy
		         max_its=atoi(optarg);
	             break;

	           case 'h': // help
	        	   printf("Usage:- -e <energy> \
	        			   -n <Number of objects> \
	        			   -g <type of gravity> 0=constant, 1=normal(inverse square), 2=inverse, 3=\"SHM\", 4=super gravity \n");
		           exit(0);
	             break;
	           }
    }

	screen_diag=sqrt(H*H+W*W);

	#ifdef SDL_ON
    /* Initialize SDL */
    SDL_Init(SDL_INIT_VIDEO);

    /* Initialize the screen / window */
    screen = SDL_SetVideoMode(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_DEPTH, SDL_SWSURFACE|SDL_FULLSCREEN);

#endif


#ifdef ALLEGRO_ON
		allegro_init();
		install_timer();
		install_mouse();
		install_keyboard();
		srand(time(NULL));
		set_gfx_mode(MODE, W, H, 0, 0);
		buffer = create_bitmap(W, H);
#endif
		int loop_count;
		do {
			if (array_size > MAX_ARRAY_SIZE){
						array_size = MAX_ARRAY_SIZE;
			}
			if (array_size < 2){
				array_size=2;
			}

			if (world_type > 6) {
					world_type=0;
			}
		    if (world_type < 0){
		    	world_type=6;
		    }
			loop_count++;
		} while (do_world(world_type,array_size));

#ifdef ALLEGRO_ON
		allegro_exit();
#endif
#ifdef SDL_ON
		SDL_Quit();
#endif
  return 0;


}

#ifdef ALLEGRO_ON
END_OF_MAIN()
#endif

