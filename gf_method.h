/*
 * gf_method.h
 *
 * Parses argv to figure out the flags and arguments.  Creates the gf.
 */

#pragma once

#include "gf.h"

/* This prints out the error string defining the methods that you can put on argv*/
extern void methods_to_stderr();

/* Parses argv starting at "starting" */
extern int create_gf_from_argv(gf_t *gf, int w, int argc, char **argv, int starting);
