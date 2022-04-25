/////////////////////////////////////////////////////////////
// ldpc.c : an interactive, naive LDPC decoder
//
// (c) 2022 Mike Field <hamster@snap.net.nz>
//
// A product/sum LDPC decoder, following the example in
// "Introducing Low-Density Parity-Check Codes' by 
// Sarah J. Johnson.
//
// This is not at all a smart or fast implementation as it
// was written to understand the algoritim. It has been a
// big help for me. Maybe it will be a big help for others.
//
/////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <ncurses.h>
#include <math.h>

// LDPC matrix
const uint8_t matrix[4][6] = {
      {1, 1, 0, 1, 0, 0}, 
      {0, 1, 1, 0, 1, 0}, 
      {1, 0, 0, 0, 1, 1}, 
      {0, 0, 1, 1, 0, 1}
};

// Initial channel values as per page 38
double initial_r[] = { -0.5, 2.50, -4.0, 5.0, -3.5, 2.5};

#define N_ITERATIONS 8

const static char *welcome_msg[] = {
  "ldpc_demo : A simple LDPC decoder",
  "",
  "by Mike Field <hamster@snap.net.nz>",
  "",
  "LDPC codes are hard to get started with, well it was for me.",
  "",
  "This is an implementation of the example found in 'Introducing",
  "Low-Density Parity-Check Codes' by Sarah J. Johnson",
  "",
  "Keys are:",
  "  Up/Down     - change the input probablilty for the current bit",
  "  Left/Right  - Select the prior or next bit for changin",
  "  PgUp/PgDown - View the different iterations of the LDPC decoder.",
  "  ESC or Q    - Quit",
  "",
  "Hope this comes in useful for somebody. If so, send me an email!",
  "",
  "Press enter to continue:"
};

struct iteration {
   struct iteration *next;
   double **message_v_to_c;
   double *l;
   double **message_c_to_v;
   uint8_t *codeword;
   uint8_t *parity;
};

struct state {
   int cursor;
   int page;
   int iterations;
   int n_v;
   int n_c;
   double *channel;
   double *channel_llr;
   uint8_t **matrix;
   struct iteration *first_iteration;
};

// Little helper functions
static double l_to_p(double l) {
   return exp(l)/(1+exp(l));
}


static double p_to_l(double p) {
   return log(p/(1-p));
}

// A very basic display function
static void state_display(struct state *s) {
   int line = 0;
   move(line, 0);
   attron(COLOR_PAIR(1));
   printw("Channel");
   attron(COLOR_PAIR(2));
   line++;
   for(int i = 0; i < s->n_v; i++) {
      move(line,i*8);
      printw("%7.4f",s->channel[i]);
   }
   line++;

   move(line, 0);
   attron(COLOR_PAIR(1));
   printw("Channel LLR");
   attron(COLOR_PAIR(2));
   line++;
   for(int i = 0; i < s->n_v; i++) {
      move(line,i*8);
      printw("%7.4f ",s->channel_llr[i]);
   }
   line++;

   struct iteration *current = s->first_iteration;
   for(int i = 0; i < s->page; i++) {
      if(current != NULL)
        current = current->next;
   }

   if(current != NULL) {
      line++;
      move(line, 0);
      attron(COLOR_PAIR(3));
      printw("Iteraton %d of %d:", s->page+1, s->iterations);
      line++;

      move(line, 0);
      attron(COLOR_PAIR(1));
      printw("Check-to-value messages:");
      attron(COLOR_PAIR(2));
      line++;
      for(int i = 0; i < s->n_c; i++) {
         for(int j = 0; j < s->n_v; j++) {
            if(s->matrix[i][j]) {
               move(line+i, j*8);
               printw("%7.4f ",current->message_c_to_v[i][j]);
            }
         }
      }
      line += s->n_c;

      line++;

      move(line, 0);
      attron(COLOR_PAIR(1));
      printw("Value-to-check messages:");
      attron(COLOR_PAIR(2));
      line++;
      for(int i = 0; i < s->n_c; i++) {
         for(int j = 0; j < s->n_v; j++) {
            if(s->matrix[i][j]) {
               move(line+i, j*8);
               printw("%7.4f ",current->message_v_to_c[i][j]);
            }
         }
      }
      line += s->n_c;

      move(line, 0);
      attron(COLOR_PAIR(1));
      printw("L:");
      attron(COLOR_PAIR(2));
      line++;
      for(int i = 0; i < s->n_v; i++) {
         move(line,i*8);
         printw("%7.4f ",current->l[i]);
      }

      line++;

      move(line, 0);
      attron(COLOR_PAIR(1));
      printw("Codeword:");
      attron(COLOR_PAIR(2));
      line++;
      for(int i = 0; i < s->n_v; i++) {
         move(line,i*2);
         printw("%c", (current->codeword[i]) ? '1' : '0');
      }
      line++;

      move(line, 0);
      attron(COLOR_PAIR(1));
      printw("Parity:");
      attron(COLOR_PAIR(2));
      line++;
      int valid = 1;
      for(int i = 0; i < s->n_c; i++) {
         move(line,i*2);
         printw("%c", (current->parity[i]) ? '1' : '0');
         if(current->parity[i])
            valid = 0;
      }
      line++;

      line++;

      move(line,1);
      attron(valid ? COLOR_PAIR(5) : COLOR_PAIR(4));
      printw("=== %s ===   ", valid ? " Valid codeword " : "Invalid codeword");
      line++;

      current = current->next;
   }
   move(1,s->cursor*8+6);
   refresh();
}


static void add_iteration(struct state *s) {
   // Assumes malloc() always succeeds...
   struct iteration *new_i = malloc(sizeof(struct iteration));
   new_i->l                = malloc(sizeof(double)  *s->n_v);
   new_i->message_v_to_c   = malloc(sizeof(double *)*s->n_c);
   new_i->message_c_to_v   = malloc(sizeof(double *)*s->n_c);
   new_i->codeword         = malloc(sizeof(uint8_t) *s->n_v);
   new_i->parity           = malloc(sizeof(uint8_t) *s->n_c);
   for(int i = 0; i < s->n_c; i++) { 
      new_i->message_v_to_c[i] = malloc(sizeof(double)*s->n_v);
      new_i->message_c_to_v[i] = malloc(sizeof(double)*s->n_v);
   }
   // Add to list
   new_i->next = s->first_iteration;
   s->first_iteration = new_i;
}


static struct state *state_new(int n_i) {
   // Assumes malloc() always succeeds...
   int n_v = sizeof(matrix[0])/sizeof(matrix[0][0]);
   int n_c = sizeof(matrix)/sizeof(matrix[0]);
   struct state *s = malloc(sizeof(struct state));
   s->channel      = malloc(sizeof(double) * n_v);
   s->channel_llr  = malloc(sizeof(double) * n_v);
   s->matrix       = malloc(sizeof(uint8_t *) * n_c);
   for(int i = 0; i < n_c; i++) {
      s->matrix[i] = malloc(sizeof(uint8_t) * n_v); 
   }

   // Setting all the elements
   s->n_v             = n_v;
   s->n_c             = n_c;
   s->cursor          = 0; 
   s->page            = 0;
   s->iterations      = n_i;
   s->first_iteration = NULL;

   // Set the initial channel probabilities
   for(int i = 0; i < n_v; i++) {
      if(i < sizeof(initial_r)/sizeof(double)) {
         s->channel[ i] = l_to_p(initial_r[i]);
      } else {
         s->channel[ i] = 0.50;
      } 
   }
   // Copy the LDPC matrix into the structure
   for(int i = 0; i < n_v; i++) {
      for(int j = 0; j < n_c; j++) {
         s->matrix[j][i] = matrix[j][i];
      }
   }

   // Add the storage needed for each iteration 
   for(int i = 0; i < n_i; i++) {
     add_iteration(s);
   }
   return s;
}


double calc_message_v_to_c(struct state *s, struct iteration *iteration, int v, int c) {
   double t = 1.0;
   for(int i = 0; i < s->n_v; i++) {
      if(i != v && s->matrix[c][i]) {
         t *= tanh(iteration->message_c_to_v[c][i]/2);
      }
   }
   return log((1+t)/(1-t));
}

double calc_message_c_to_v(struct state *s, struct iteration *iteration, int v, int c) {
   double p = s->channel_llr[v];
   for(int i = 0; i < s->n_c; i++) {
      if(i != c && s->matrix[i][v] == 1) {
         p += iteration->message_v_to_c[i][v];
      }
   }
   return p;
} 


void state_solve(struct state *s) {
   for(int i = 0; i < s->n_v; i++) {
      s->channel_llr[i] = p_to_l(s->channel[i]);
   }
   if(s->first_iteration == NULL) 
     return;

   for(int c = 0; c < s->n_c; c++) {  // For each check node
      for(int v = 0; v < s->n_v; v++) {  // message for each check node 
         s->first_iteration->message_c_to_v[c][v] = s->channel_llr[v];
      }
   }

   struct iteration *current = s->first_iteration;
   while(current != NULL) {
      for(int c = 0; c < s->n_c; c++) {  // For each check node
         for(int v = 0; v < s->n_v; v++) {  // message for each check node 
            current->message_v_to_c[c][v] = calc_message_v_to_c(s, current, v, c);
         }
      }

      for(int i = 0; i < s->n_c; i++) {
         current->parity[i] = 0;
      }
      for(int v = 0; v < s->n_v; v++) {
         int b;
         current->l[v] = s->channel_llr[v];
         for(int i = 0; i < s->n_c; i++) {
            if(s->matrix[i][v])
              current->l[v] += current->message_v_to_c[i][v];
         }
         b = (current->l[v] < 0) ? 1 : 0;
         current->codeword[v] = b;
         for(int i = 0; i < s->n_c; i++) {
            if(s->matrix[i][v])
               current->parity[i] ^= b;
         }
         
      }
      
      if(current->next != NULL) {
         for(int c = 0; c < s->n_c; c++) {  // For each check node
            for(int v = 0; v < s->n_v; v++) {  // message for each check node 
                current->next->message_c_to_v[c][v] = calc_message_c_to_v(s, current, v, c);
            }
         }
      }
      current = current->next;
   }
}

void state_delete(struct state *s) {
   while(s->first_iteration != NULL) {
     struct iteration *x = s->first_iteration;
     s->first_iteration = x->next;
     free(x->l);
     free(x->codeword);
     free(x->parity);
     for(int i = 0; i < s->n_c; i++) {
        free(x->message_v_to_c[i]);
        free(x->message_c_to_v[i]);
     }
     free(x->message_v_to_c);
     free(x->message_c_to_v);
     free(x);
   }
   free(s->channel);
   free(s->channel_llr);
   for(int i = 0; i < s->n_c; i++) {
      free(s->matrix[i]); 
   }
   free(s->matrix);
   free(s);
}

int process_keys(struct state *s) {
   int key = getch();
   switch(key) {
      case 27  :  return 0; 

      case 'q' :
      case 'Q' :  return 0; 

      case KEY_PPAGE:  if(s->page > 0)
                          s->page--;
                       break;
 
      case KEY_NPAGE:  if(s->page < s->iterations-1)
                          s->page++;
                       break;
 
      case KEY_LEFT:   if(s->cursor == 0)
                           s->cursor = s->n_v-1;
                       else
                           s->cursor--;
                       break;
 
      case KEY_RIGHT:  if(s->cursor == s->n_v-1)
                           s->cursor = 0;
                       else
                           s->cursor++;
                       break;

      case KEY_UP:     s->channel[s->cursor] = floor(s->channel[s->cursor]*100)/100; 
                       if(s->channel[s->cursor] > 0.01)
                          s->channel[s->cursor] -= 0.01;
                       else
                          s->channel[s->cursor] = 0.01;
                       state_solve(s);
                       break;

      case KEY_DOWN:   s->channel[s->cursor] = floor(s->channel[s->cursor]*100)/100;
                       if(s->channel[s->cursor] < 0.99)
                          s->channel[s->cursor] += 0.01;
                       else
                          s->channel[s->cursor] = 0.99;
                       state_solve(s);

                       break;
   };

   return 1;
}

void welcome_screen(void) {
  clear();
  attron(COLOR_PAIR(2));
  for(int i = 0; i < sizeof(welcome_msg)/sizeof(char*); i++) {
     move(i,0);
     printw(welcome_msg[i]);
     attron(COLOR_PAIR(1));
  }
  getch();
  clear();
}


int main(int argc, char *argv[]) {
   struct state *s;
   initscr();
   if(!has_colors()) {
      endwin();
      fprintf(stderr,"Console does not support color\n");
      return 0;
   }
   start_color();
   init_pair(1, COLOR_GREEN, COLOR_BLACK);  // Headings
   init_pair(2, COLOR_WHITE, COLOR_BLACK);  // General text
   init_pair(3, COLOR_RED,   COLOR_BLACK);  // Section heading
   init_pair(4, COLOR_WHITE, COLOR_RED);    // Invalid codeword message
   init_pair(5, COLOR_WHITE, COLOR_GREEN);  // Valid codeword message

   keypad(stdscr,TRUE);

   welcome_screen(); 

   s = state_new(N_ITERATIONS);
   state_solve(s);

   do {
      state_display(s);
   } while(process_keys(s));
   state_delete(s);

   endwin();
}
