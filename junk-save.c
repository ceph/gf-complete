
/*
      c = gf.multiply.w32(&gf, a, b);
      tested = 0;

*/
      /* If this is not composite, then first test against the default: */

/*
      if (h->mult_type != GF_MULT_COMPOSITE) {
        tested = 1;
        d = gf_def.multiply.w32(&gf_def, a, b);

        if (c != d) {
          printf("Error in single multiplication (all numbers in hex):\n\n");
          printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
          printf("  The default returned %x\n", d);
          exit(1);
        }
      }

*/
      /* Now, we also need to double-check, in case the default is wanky, and when
         we're performing composite operations. Start with 0 and 1: */

/*
      if (a == 0 || b == 0 || a == 1 || b == 1) {
        tested = 1;
        if (((a == 0 || b == 0) && c != 0) ||
            (a == 1 && c != b) ||
            (b == 1 && c != a)) {
          printf("Error in single multiplication (all numbers in hex):\n\n");
          printf("  gf.multiply.w32() of %x and %x returned %x, which is clearly wrong.\n", a, b, c);
          exit(1);
        }
      
*/
      /* If division or inverses are defined, let's test all combinations to make sure
         that the operations are consistent with each other. */
          
/*
      } else {
        if ((c & mask) != c) {
          printf("Error in single multiplication (all numbers in hex):\n\n");
          printf("  gf.multiply.w32() of %x and %x returned %x, which is too big.\n", a, b, c);
          exit(1);
        }
          
      }
      if (gf.inverse.w32 != NULL && (a != 0 || b != 0)) {
        tested = 1;
        if (a != 0) {
          ai = gf.inverse.w32(&gf, a);

          if (gf.multiply.w32(&gf, c, ai) != b) {
            printf("Error in single multiplication (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
            printf("  The inverse of %x is %x, and gf_multiply.w32() of %x and %x equals %x\n",
                      a, ai, c, ai, gf.multiply.w32(&gf, c, ai));
            exit(1);
          }
        }
        if (b != 0) {
          bi = gf.inverse.w32(&gf, b);
          if (gf.multiply.w32(&gf, c, bi) != a) {
            printf("Error in single multiplication (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
            printf("  The inverse of %x is %x, and gf_multiply.w32() of %x and %x equals %x\n",
                      b, bi, c, bi, gf.multiply.w32(&gf, c, bi));
            exit(1);
          }
        }
      }
      if (gf.divide.w32 != NULL && (a != 0 || b != 0)) {
        tested = 1;
      
        if (a != 0) {
          ai = gf.divide.w32(&gf, c, a);

          if (ai != b) {
            printf("Error in single multiplication/division (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
            printf("  gf.divide.w32() of %x and %x returned %x\n", c, a, ai);
            exit(1);
          }
        }
        if (b != 0) {
          bi = gf.divide.w32(&gf, c, b);

          if (bi != a) {
            printf("Error in single multiplication/division (all numbers in hex):\n\n");
            printf("  gf.multiply.w32() of %x and %x returned %x\n", a, b, c);
            printf("  gf.divide.w32() of %x and %x returned %x\n", c, b, bi);
            exit(1);
          }
        }
      }

      if (!tested) problem("There is no way to test multiplication.\n");
    }
*/

/*
  if (region) {
    
    if (w == 4) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r8b = (uint8_t *) malloc(REGION_SIZE);
        r8c = (uint8_t *) malloc(REGION_SIZE);
        r8d = (uint8_t *) malloc(REGION_SIZE);
        fill_random_region(r8b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 16; a++) {
            fill_random_region(r8c, REGION_SIZE);
            memcpy(r8d, r8c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8c+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = (r8b[i] >> 4);
              c = (r8c[i] >> 4);
              d = (r8d[i] >> 4);
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d * %d = %d, but should equal %d\n", a, b, c, gf.multiply.w32(&gf, a, b) );
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d %d %d %d\n", a, b, c, d);
                printf("       %d %d %d %d\n", a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
              b = (r8b[i] & 0xf);
              c = (r8c[i] & 0xf);
              d = (r8d[i] & 0xf);
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       %d * %d = %d, but should equal %d\n", a, b, c, gf.multiply.w32(&gf, a, b) );
                printf("i=%d.  0x%x 0x%x 0x%x 0x%x\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r8b+i));
                printf("       (%d * %d ^ %d) should equal %d - equals %d\n", 
                    a, b, d, (gf.multiply.w32(&gf, a, b) ^ d), c);
                printf("       %d %d %d %d\n", a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 16; a++) {
            fill_random_region(r8b, REGION_SIZE);
            memcpy(r8d, r8b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8b+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = (r8b[i] >> 4);
              d = (r8d[i] >> 4);
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("i=%d.  %d %d %d\n", i, a, b, d);
                printf("i=%d.  %d %d %d\n", i, a, r8b[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
              b = (r8b[i] & 0xf);
              d = (r8d[i] & 0xf);
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("%d %d %d\n", a, b, d);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r8b);
        free(r8c);
        free(r8d);
      }
    } else if (w == 8) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r8b = (uint8_t *) malloc(REGION_SIZE);
        r8c = (uint8_t *) malloc(REGION_SIZE);
        r8d = (uint8_t *) malloc(REGION_SIZE);
        fill_random_region(r8b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 256; a++) {
            fill_random_region(r8c, REGION_SIZE);
            memcpy(r8d, r8c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE;
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8c+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w8(i, (uint8_t*)r8b, REGION_SIZE / 2);
                c = get_alt_map_2w8(i, (uint8_t*)r8c, REGION_SIZE / 2);
                d = get_alt_map_2w8(i, (uint8_t*)r8d, REGION_SIZE / 2);
              } else {
                b = r8b[i];
                c = r8c[i];
                d = r8d[i];
              }
              if (!xor && gf.multiply.w32(&gf, a, b) != c) {
                printf("i=%d.  %d %d %d %d\n", i, a, b, c, d);
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                printf("%llx.  Sindex: %d\n", r8b+i, sindex);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w32(&gf, a, b) ^ d) != c) {
                printf("i=%d.  %d %d %d %d\n", i, a, b, c, d);
                printf("i=%d.  %d %d %d %d\n", i, a, r8b[i], r8c[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;  
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (a = 0; a < 256; a++) {
            fill_random_region(r8b, REGION_SIZE);
            memcpy(r8d, r8b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint8_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint8_t);
            gf.multiply_region.w32(&gf, (void *) (r8b+sindex), (void *) (r8b+sindex), a, size, xor);
            for (i = sindex; i < eindex; i++) {
              b = r8b[i];
              d = r8d[i];
              if (!xor && gf.multiply.w32(&gf, a, d) != b) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w32(&gf, a, d) ^ d) != b) {
                printf("i=%d.  %d %d %d\n", i, a, b, d);
                printf("i=%d.  %d %d %d\n", i, a, r8b[i], r8d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r8b);
        free(r8c);
        free(r8d);
      }
    } else if (w == 16) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r16b = (uint16_t *) malloc(REGION_SIZE);
        r16c = (uint16_t *) malloc(REGION_SIZE);
        r16d = (uint16_t *) malloc(REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            fill_random_region(r16b, REGION_SIZE);
            a = MOA_Random_W(w, 0);
            fill_random_region(r16c, REGION_SIZE);
            memcpy(r16d, r16c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE / sizeof(uint16_t);
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint16_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint16_t);
            gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16c+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r16c+sindex), (void *) (r16d+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r16c+sindex), (void *) (r16d+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r16d+sindex), (void *) (r16b+sindex), ai, size, xor);
            }
  
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w16(i, (uint8_t*)r16b, size / 2);
                c = get_alt_map_2w16(i, (uint8_t*)r16c, size / 2);
                d = get_alt_map_2w16(i, (uint8_t*)r16d, size / 2);
              } else {
                b = r16b[i];
                c = r16c[i];
                d = r16d[i];
              }
              if (!xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We have %d * %d = %d, and %d * %d = %d.\n", a, b, c, c, ai, d);
                printf("%d is the inverse of %d\n", ai, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && b != 0) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("   b should equal 0, but it doesn't.  Probe into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_W(w, 0);
            fill_random_region(r16b, REGION_SIZE);
            memcpy(r16d, r16b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint16_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint16_t);
            gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r16d+sindex), (void *) (r16b+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r16b+sindex), (void *) (r16b+sindex), ai, size, 0);
            }
  
            for (i = sindex; i < eindex; i++) {
              b = r16b[i];
              c = r16c[i];
              d = r16d[i];
              if (!xor && (d != b)) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=b; b = ba; b = b(a^-1).\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                printf("b = %d.  d = %d.  a = %d\n", b, d, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r16b+i));
                printf("We did d=b; b = b + ba; b += d; b = b(a^-1);\n");
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r16b);
        free(r16c);
        free(r16d);
      }
    } else if (w == 32) {
      if (gf.multiply_region.w32 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r32b = (uint32_t *) malloc(REGION_SIZE);
        r32c = (uint32_t *) malloc(REGION_SIZE);
        r32d = (uint32_t *) malloc(REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_32();
            fill_random_region(r32b, REGION_SIZE);
            fill_random_region(r32c, REGION_SIZE);
            memcpy(r32d, r32c, REGION_SIZE);
            if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
              sindex = 0;
              eindex = REGION_SIZE / sizeof(uint32_t);
            } else {
              sindex = MOA_Random_W(3, 1);
              eindex = REGION_SIZE/sizeof(uint32_t)-MOA_Random_W(3, 1);
            }
            size = (eindex-sindex)*sizeof(uint32_t);
            gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32c+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r32c+sindex), (void *) (r32d+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r32c+sindex), (void *) (r32d+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r32d+sindex), (void *) (r32b+sindex), ai, size, xor);
            }
            for (i = sindex; i < eindex; i++) {
              if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
                (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
                b = get_alt_map_2w32(i, (uint8_t*)r32b, size / 2);
                c = get_alt_map_2w32(i, (uint8_t*)r32c, size / 2);
                d = get_alt_map_2w32(i, (uint8_t*)r32d, size / 2);
                i++;
              } else {
                b = r32b[i];
                c = r32c[i];
                d = r32d[i];
              }
              if (!xor && d != b) {
                printf("i=%d.  Addresses: b: 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We have %d * %d = %d, and %d * %d = %d.\n", a, b, c, c, ai, d);
                printf("%d is the inverse of %d\n", ai, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && b != 0) {
                printf("i=%d.  Addresses: b: 0x%lx   c: 0x%lx   d: 0x%lx\n", i, 
                      (unsigned long) (r32b+i), (unsigned long) (r32c+i), (unsigned long) (r32d+i));
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("   b should equal 0, but it doesn't.  Probe into it.\n");
                printf("a: %8x  b: %8x  c: %8x,  d: %8x\n", a, b, c, d);
                problem("Failed buffer-constant, xor=1");
              }

            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if ((((gf_internal_t*)gf.scratch)->region_type & GF_REGION_ALTMAP) != 0 &&
            (((gf_internal_t*)gf.scratch)->mult_type == GF_MULT_COMPOSITE)) {
            continue;
          }
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a = MOA_Random_32();
            fill_random_region(r32b, REGION_SIZE);
            memcpy(r32d, r32b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint32_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint32_t);
            gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), a, size, xor);
            ai = gf.inverse.w32(&gf, a);
            if (!xor) {
              gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), ai, size, xor);
            } else {
              gf.multiply_region.w32(&gf, (void *) (r32d+sindex), (void *) (r32b+sindex), 1, size, xor);
              gf.multiply_region.w32(&gf, (void *) (r32b+sindex), (void *) (r32b+sindex), ai, size, 0);
            }

            for (i = sindex; i < eindex; i++) {
              b = r32b[i];
              c = r32c[i];
              d = r32d[i];
              if (!xor && (d != b)) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We did d=b; b = ba; b = b(a^-1).\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                printf("b = %d.  d = %d.  a = %d\n", b, d, a);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && d != b) {
                printf("i=%d.  Address 0x%lx\n", i, (unsigned long) (r32b+i));
                printf("We did d=b; b = b + ba; b += d; b = b(a^-1);\n");
                printf("We did d=c; c ^= ba; d ^= c; b ^= (a^-1)d;\n");
                printf("So, b should equal d, but it doesn't.  Look into it.\n");
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r32b);
        free(r32c);
        free(r32d);
      }
    } else if (w == 64) {
      if (gf.multiply_region.w64 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r64b = (uint64_t *) malloc(REGION_SIZE);
        r64c = (uint64_t *) malloc(REGION_SIZE);
        r64d = (uint64_t *) malloc(REGION_SIZE);
        fill_random_region(r64b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a64 = MOA_Random_64();
            fill_random_region(r64c, REGION_SIZE);
            memcpy(r64d, r64c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint64_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t);
            gf.multiply_region.w64(&gf, (void *) (r64b+sindex), (void *) (r64c+sindex), a64, size, xor);
            for (i = sindex; i < eindex; i++) {
              b64 = r64b[i];
              c64 = r64c[i];
              d64 = r64d[i];
              if (!xor && gf.multiply.w64(&gf, a64, b64) != c64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx should be 0x%llx\n", i, a64, b64, c64, 
                            gf.multiply.w64(&gf, a64, b64));
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64c[i]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && (gf.multiply.w64(&gf, a64, b64) ^ d64) != c64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx 0x%llx\n", i, a64, b64, c64, d64);
                printf("i=%d.  0x%llx 0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64c[i], r64d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            a64 = MOA_Random_64();
            fill_random_region(r64b, REGION_SIZE);
            memcpy(r64d, r64b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/sizeof(uint64_t)-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t);
            gf.multiply_region.w64(&gf, (void *) (r64b+sindex), (void *) (r64b+sindex), a64, size, xor);
            for (i = sindex; i < eindex; i++) {
              b64 = r64b[i];
              d64 = r64d[i];
              if (!xor && gf.multiply.w64(&gf, a64, d64) != b64) problem("Failed buffer-constant, xor=0");
              if (xor && (gf.multiply.w64(&gf, a64, d64) ^ d64) != b64) {
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, b64, d64);
                printf("i=%d.  0x%llx 0x%llx 0x%llx\n", i, a64, r64b[i], r64d[i]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r64b);
        free(r64c);
        free(r64d);
      }
    } else if (w == 128) {
      if (gf.multiply_region.w128 == NULL) {
        printf("No multiply_region.\n");
      } else {
        r128b = (uint64_t *) malloc(REGION_SIZE);
        r128c = (uint64_t *) malloc(REGION_SIZE);
        r128d = (uint64_t *) malloc(REGION_SIZE);
        fill_random_region(r128b, REGION_SIZE);
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src != dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            MOA_Random_128(a128);
            fill_random_region(r128c, REGION_SIZE);
            memcpy(r128d, r128c, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            eindex = REGION_SIZE/(2*sizeof(uint64_t))-MOA_Random_W(3, 1);
            size = (eindex-sindex)*sizeof(uint64_t)*2;
            gf.multiply_region.w128(&gf, (void *) (r128b+sindex*2), (void *) (r128c+sindex*2), a128, size, xor);
            for (i = sindex; i < eindex; i++) {
              b128[0] = r128b[2*i];
              b128[1] = r128b[2*i+1];
              c128[0] = r128c[2*i];
              c128[1] = r128c[2*i+1];
              d128[0] = r128d[2*i];
              d128[1] = r128d[2*i+1];
              gf.multiply.w128(&gf, a128, b128, e128);
              if (xor) {
                e128[0] ^= d128[0];
                e128[1] ^= d128[1];
              }
              if (!xor && !GF_W128_EQUAL(c128, e128)) {
                printf("i=%d.  0x%llx%llx 0x%llx%llx 0x%llx%llx should be 0x%llx%llx\n",
                        i, a128[0], a128[1], b128[0], b128[1], c128[0], c128[1], e128[0], e128[1]);
                problem("Failed buffer-constant, xor=0");
              }
              if (xor && !GF_W128_EQUAL(e128, c128)) {
                printf("i=%d.  0x%llx%llx 0x%llx%llx 0x%llx%llx 0x%llx%llx\n", i,
                        a128[0], a128[1], b128[0], b128[1], c128[0], c128[1], d128[0], d128[1]);
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        for (xor = 0; xor < 2; xor++) {
          if (verbose) {
            printf("Testing buffer-constant, src == dest, xor = %d\n", xor);
            fflush(stdout);
          }
          for (j = 0; j < 1000; j++) {
            MOA_Random_128(a128);
            fill_random_region(r128b, REGION_SIZE);
            memcpy(r128d, r128b, REGION_SIZE);
            sindex = MOA_Random_W(3, 1);
            sindex = 0;
            eindex = REGION_SIZE/(2*sizeof(uint64_t))-MOA_Random_W(3, 1);
            eindex = REGION_SIZE/(2*sizeof(uint64_t));
            size = (eindex-sindex)*sizeof(uint64_t)*2;
            gf.multiply_region.w128(&gf, (void *) (r128b+sindex), (void *) (r128b+sindex), a128, size, xor);
            for (i = sindex; i < eindex; i++) {
              b128[0] = r128b[2*i];
              b128[1] = r128b[2*i + 1];
              d128[0] = r128d[2*i];
              d128[1] = r128d[2*i + 1];
              gf.multiply.w128(&gf, a128, d128, e128);
              if (xor) {
                e128[0] ^= d128[0];
                e128[1] ^= d128[1];
              }
              if (!xor && !GF_W128_EQUAL(b128, e128)) problem("Failed buffer-constant, xor=0");
              if (xor && !GF_W128_EQUAL(b128, e128)) {
                problem("Failed buffer-constant, xor=1");
              }
            }
          }
        }
        free(r128b);
        free(r128c);
        free(r128d);
      }
    }
  }
  exit(0);
  */
}
