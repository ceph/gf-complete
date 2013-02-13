/*
 * gf_poly.c - program to help find primitive polynomials in composite fields
 */

#include "gf_complete.h"
#include "gf_method.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GF_POLY_COEF_MASK8 0xff
#define GF_POLY_COEF_MASK16 0xffff
#define GF_POLY_COEF_MASK32 0xffffffff
#define GF_POLY_COEF_MASK64 0xffffffffffffffff

#define LLUI (long long unsigned int)

struct gf_poly_coef_s;

typedef struct gf_poly_coef_s {
  uint64_t coef;
  uint64_t power;
  struct gf_poly_coef_s *next;
} gf_poly_coef_t;

typedef struct gf_poly_s {
  gf_poly_coef_t *leading_coef;
  uint64_t num_coefs;
  gf_t *coef_gf;
  int w;
} gf_poly_t;

static uint64_t gf_add(int w, uint64_t a, uint64_t b)
{
  if (w == 8) {
    return (a & GF_POLY_COEF_MASK8) ^ (b & GF_POLY_COEF_MASK8);
  } else if (w == 16) {
    return (a & GF_POLY_COEF_MASK16) ^ (b & GF_POLY_COEF_MASK16);
  } else if (w == 32) {
    return (a & GF_POLY_COEF_MASK32) ^ (b & GF_POLY_COEF_MASK32);
  } else if (w == 64) {
    return (a & GF_POLY_COEF_MASK64) ^ (b & GF_POLY_COEF_MASK64);
  }
}

static uint64_t gf_mult(int w, gf_t* gf, uint64_t a, uint64_t b)
{
  if (w <= 32) {
    return gf->multiply.w32(gf, a, b); 
  } else if (w == 64) {
    return gf->multiply.w64(gf, a, b); 
  }
}

static uint64_t gf_divide(int w, gf_t* gf, uint64_t a, uint64_t b)
{
  if (w <= 32) {
    return gf->divide.w32(gf, a, b); 
  } else if (w == 64) {
    return gf->divide.w64(gf, a, b); 
  }
}

static uint64_t gf_inverse(int w, gf_t* gf, uint64_t a)
{
  if (w <= 32) {
    return gf->inverse.w32(gf, a); 
  } else if (w == 64) {
    return gf->inverse.w64(gf, a);
  }
}

gf_poly_t* gf_poly_init(int w, gf_t *gf)
{
  gf_poly_t *gf_poly = (gf_poly_t*)malloc(sizeof(gf_poly_t));

  if (gf_poly == NULL || gf == NULL) {
    return NULL;
  }

  gf_poly->leading_coef = NULL;
  gf_poly->num_coefs = 0;
  gf_poly->coef_gf = gf;
  gf_poly->w = w;

  return gf_poly;
}

void gf_poly_print(gf_poly_t *gf_poly, char *message)
{
  gf_poly_coef_t *tmp;

  if (gf_poly == NULL) {
    fprintf(stderr, "0 * x^0\n");
    return;
  }

  tmp = gf_poly->leading_coef;

  while (tmp != NULL) {
    printf("%llu * x^%llu", LLUI tmp->coef, LLUI tmp->power);
    tmp = tmp->next;
    if (tmp) {
      printf(" + ");
    }
  }

  if (message != NULL) {
    printf(": %s\n", message);
  }
}

gf_poly_t* gf_poly_copy(gf_poly_t *poly)
{
  gf_poly_t *new_poly = (gf_poly_t*)malloc(sizeof(gf_poly_t));
  gf_poly_coef_t *tmp = poly->leading_coef;

  if (new_poly == NULL) {
    return NULL;
  }

  new_poly->leading_coef = NULL;
  new_poly->num_coefs = 0;
  new_poly->coef_gf = poly->coef_gf;
  new_poly->w = poly->w;
  
  while (tmp != NULL) {
    gf_poly_add_coef(new_poly, tmp->coef, tmp->power);

    tmp = tmp->next;
  }

  return new_poly;
}

void gf_poly_clear(gf_poly_t* a)
{
  while (a->leading_coef != NULL) {
    gf_poly_coef_t *tmp = a->leading_coef;
    
    a->leading_coef = tmp->next;

    free(tmp);
  }
}

void gf_poly_free(gf_poly_t **a)
{
  gf_poly_clear(*a);
  free(*a); 
  *a = NULL;
}

gf_poly_coef_t* gf_poly_create_node(uint64_t coef, uint64_t power)
{
  gf_poly_coef_t* node = (gf_poly_coef_t*)malloc(sizeof(gf_poly_coef_t));

  if (node == NULL) {
    return NULL;
  }

  node->coef = coef;
  node->power = power;
  node->next = NULL;

  return node;
}

int gf_poly_remove_node(gf_poly_t *gf_poly, uint64_t power)
{
  gf_poly_coef_t* iter = gf_poly->leading_coef;

  if (iter->power == power) {
    gf_poly->leading_coef = iter->next;   
    free(iter);
    return 0;
  }

  while (iter->next != NULL) {
    if (iter->next->power == power) {
      gf_poly_coef_t* tmp = iter->next;
      iter->next = iter->next->next;
      free(tmp);
      return 0;
    }
    iter = iter->next;
  }

  return -1;
}

int gf_poly_add_coef(gf_poly_t *gf_poly, uint64_t coef_val, uint64_t power)
{
  gf_poly_coef_t* node;
  gf_poly_coef_t* iter = gf_poly->leading_coef;

  /*
   * The new node has the highest power, or there are no terms
   */
  if (gf_poly->leading_coef == NULL || gf_poly->leading_coef->power < power) {
    node = gf_poly_create_node(coef_val, power);
    node->next = gf_poly->leading_coef;
    gf_poly->leading_coef = node;
    return 0;
  }

  /*
   * The new node is of the same power, add the coefs
   */
  if (gf_poly->leading_coef->power == power) {
    gf_poly->leading_coef->coef = gf_add(gf_poly->w, gf_poly->leading_coef->coef, coef_val);   
    if (gf_poly->leading_coef->coef == 0) {
      gf_poly_remove_node(gf_poly, power);
    }
    return 0;
  }

  while (iter->next != NULL) {
    if (iter->next->power == power) {
      iter->next->coef = gf_add(gf_poly->w, iter->next->coef, coef_val);   

      if (iter->next->coef == 0) {
        gf_poly_remove_node(gf_poly, power);
      }

      return 0;
    }
    if (iter->next->power < power) {
      node = gf_poly_create_node(coef_val, power);
      node->next = iter->next;
      iter->next = node;
      return 0;
    }
    iter = iter->next;
  }
  
  /*
   * The power passed in is lower than any in the existing poly
   */
  node = gf_poly_create_node(coef_val, power);
  iter->next = node;

  return 0;
}

/*
 * Compute a+b and store in a
 */
int gf_poly_add(gf_poly_t* a, gf_poly_t* b)
{
  gf_poly_coef_t* iter = b->leading_coef;

  while (iter != NULL) {
    gf_poly_add_coef(a, iter->coef, iter->power);
    iter = iter->next; 
  }

  return 0;
}

/*
 * Compute a*b and store in a
 */
int gf_poly_mult(gf_poly_t* a, gf_poly_t* b)
{
  gf_poly_coef_t* a_iter = a->leading_coef;

  /*
   * Remove one node at a time from 'a', starting with
   * highest power.  Multiply the removed (coef,power)
   * by every entry of 'b,' adding each product into 'a.'
   */
  while (a_iter != NULL) {
    gf_poly_coef_t* tmp = a_iter;
    gf_poly_coef_t* b_iter = b->leading_coef;

    uint64_t a_power = a_iter->power;
    uint64_t a_coef = a_iter->coef;
    a_iter = a_iter->next;
    gf_poly_remove_node(a, tmp->power);

    while (b_iter != NULL) {
      uint64_t new_power = b_iter->power + a_power;
      uint64_t new_coef = gf_mult(a->w, a->coef_gf, b_iter->coef, a_coef);

      gf_poly_add_coef(a, new_coef, new_power);

      b_iter = b_iter->next;
    }
  }
  return 0;
}

/*
 * Compute a % b and store in a
 */
int gf_poly_reduce(gf_poly_t* a, gf_poly_t* b)
{
   gf_poly_t* c = gf_poly_init(a->w, a->coef_gf);
   gf_poly_coef_t* a_iter = a->leading_coef;
   gf_poly_coef_t* b_iter = b->leading_coef;

  /*
   * Reduce until the degree of 'a' is less than
   * the degree of 'b.'  At that point 'a' will 
   * contain the remainder of a / b.
   */
  while (a_iter && (a_iter->power >= b_iter->power)) {

    /*
     * Get the degree and leading coef of the current
     * 'b'.
     */
    uint64_t reduce_power = a_iter->power - b_iter->power;
    uint64_t reduce_coef = gf_divide(a->w, a->coef_gf, a_iter->coef, b_iter->coef);

    /*
     * Create a poly that will get rid of leading power
     * of 'b' when added: c*x^(n-m)*b(x), where c 
     * is the leading coef of 'a', n is the deg of 'a'
     * and m is the degree of 'b'.
     */
    gf_poly_add_coef(c, reduce_coef, reduce_power);
    gf_poly_mult(c, b);
    
    /*
     * Add the newly created poly, which will reduce 
     * a(x) by at least one term (leading term).
     */
    gf_poly_add(a, c);
    
    gf_poly_clear(c); 
   
    /*
     * Grab the new leading term of 'a'
     */ 
    a_iter = a->leading_coef;
  }
}

/*
 * Get the GCD of a and b, return the result
 */
gf_poly_t* gf_poly_gcd(gf_poly_t* a, gf_poly_t* b)
{
  gf_poly_t *r1, *r2;
  gf_poly_t* tmp_swp;

  if (a->leading_coef == NULL || b->leading_coef == NULL) {
    return NULL;
  }

  if (a->leading_coef->power > b->leading_coef->power) {
    r1 = a;
    r2 = b;
  } else {
    r1 = b;
    r2 = a;
  }

  while ( 1 ) {
    if (r2->leading_coef == NULL) {
      break;
    }
    if (r2->leading_coef->power == 0 && r2->leading_coef->coef <= 1) {
      break;
    }

    gf_poly_reduce(r1, r2);
    tmp_swp = r1;
    r1 = r2;
    r2 = tmp_swp;
  }

  return r1;
}

/*
 * The Ben-Or algorithm for determining irreducibility
 */
int gf_poly_is_irred(gf_poly_t* poly)
{
  gf_poly_t *gcd;
  gf_poly_t *prod_of_irred;
  uint64_t prod_of_irred_power = ((unsigned long long) 1) << poly->w;
  int n = poly->leading_coef->power / 2;
  int i;
  int ret = 0;
  gf_poly_t *a = gf_poly_copy(poly);

  prod_of_irred = gf_poly_init(a->w, a->coef_gf);


  for (i = 1; i <= n; i++) {
    gf_poly_add_coef(prod_of_irred, 1, prod_of_irred_power);
    gf_poly_add_coef(prod_of_irred, 1, 1);
  
    gf_poly_reduce(prod_of_irred, a); 
    
    gcd = gf_poly_gcd(a, prod_of_irred); 

    /*
     * It is irreducible if it is not the product of 
     * non-trivial factors (non-constant).  Therefore,
     * the GCD of the poly and prod_of_irred should be
     * a constant (0 or 0-degree polynomial).
     */ 
    if (gcd == NULL) {
      ret = -1;
      break;
    } else if (gcd->leading_coef->power != 0) {
      ret = -1;
      break;
    } else if (gcd->leading_coef->power == 0) {
      ret = 0;
      break;
    } else {
      ret = -1;
      break;
    }
    
    // Need if to avoid a overflow error
    if ((i + 1) <= n) {
      prod_of_irred_power *= prod_of_irred_power;
    }
    gf_poly_clear(prod_of_irred);
  }

  gf_poly_free(&a);

  return ret;
}

int is_suitible_s(int w, gf_t *gf, uint64_t s)
{
  uint64_t num_elems = ((unsigned long long) 1) << w;
  uint64_t i = 2;
  uint64_t i_inv;

  for (; i < num_elems; i++) {
    i_inv = gf_inverse(w, gf, i);
    if ((i ^ i_inv) == s) {
      fprintf(stderr, "Bailed on %llu ^ %llu = %llu\n", LLUI i, LLUI i_inv, LLUI s);
      return -1;
    }
    if (i % 1000000000 == 0) fprintf(stderr, "Processed %llu\n", LLUI i);
  }

  return 0;
}

static void
usage(char *cmd)
{
  fprintf(stderr, "%s w <GF args> S <s value>\n", cmd);
  fprintf(stderr, "\t will build a trinomial x^2+S*x+1\n");
  fprintf(stderr, "OR\n");
  fprintf(stderr, "%s w <GF args> G coef1,power1 <coef2,power2> ... <coefn,powern>\n", cmd);
  fprintf(stderr, "\t will build a polynomial coef1^(power1) + ... + coefn^(powern)\n");
  fprintf(stderr, "Example: ./gf_poly 8 - - - G 1,2 2,1 1,0\n");
  fprintf(stderr, "\t will build a polynomial x^2+2*x+1 with coefs from GF(2^8)\n");
}

/*
 * Find irred poly of form x^2+sx+1
 * a_n*x^n + a_(n-1)*x^(n-1) + ...
 *
 * Terms are specified as: a_i,i a_j,j, ... where 
 * i is the degree of the term and a_i is the coef
 *
 */
int main(int argc, char **argv)
{
  gf_t gf;
  int ret;
  int w;
  int i;
  uint64_t irred_coef_s;
  gf_poly_t *irred_poly;
  char *term;

  bzero(&gf, sizeof(gf_t)); 

  if (argc < 4) {
    usage(argv[0]);
    return -1;
  }
  
  w = atoi(argv[1]);
  
  ret = create_gf_from_argv(&gf, w, argc, argv, 3);

  if (ret <= 0) {
    fprintf(stderr, "Could not create a GF\n");
    return -1;
  }
    
  irred_poly = gf_poly_init(w, &gf);

  i = ret + 1;

  if (strlen(argv[i]) > 1) {
    usage(argv[0]); 
    exit(1);
  }

  if (argv[i][0] == 'S') {
    i++;
    irred_coef_s = (uint64_t)strtoull(argv[i], NULL, 10);
  
    /*
     * If this is a trinomial of the form x^2+s*x+1, then
     * we can do a quick pre-check to see if this may be
     * an irreducible polynomial.
     */
    if (is_suitible_s(w, &gf, irred_coef_s) < 0) {
      fprintf(stderr, "%llu is not a suitable coeffient!\n", LLUI irred_coef_s);
      return -1;
    } else {
      fprintf(stderr, "%llu IS A suitable coeffient!\n", LLUI irred_coef_s);
    }


    gf_poly_add_coef(irred_poly, 1, 2);
    gf_poly_add_coef(irred_poly, irred_coef_s, 1);
    gf_poly_add_coef(irred_poly, 1, 0);

  } else if (argv[i][0] == 'G') {
    term = argv[++i];


    while (term != NULL) {
      uint64_t coef = strtoull(strtok(term, ","), NULL, 10);
      uint64_t power = strtoull(strtok(NULL, ","), NULL, 10);
    
      gf_poly_add_coef(irred_poly, coef, power);
    
      if (i < argc) {
        term = argv[++i];
      } else {
        break;
      }
    }
  } else {
    usage(argv[0]);
    exit(1);
  }
  
  gf_poly_print(irred_poly, " specified via the command line\n");

  ret = gf_poly_is_irred(irred_poly);

  if (ret < 0) {
    gf_poly_print(irred_poly, " IS NOT irreducible\n");
  } else {
    gf_poly_print(irred_poly, " IS irreducible\n");
  }

  return 0;
}
