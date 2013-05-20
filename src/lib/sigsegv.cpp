/**
 * This file is based on the sigsegv.c utility by Jaco Kroon
 * (http://tlug.up.ac.za/wiki/index.php/Obtaining_a_stack_trace_in_C_upon_SIGSEGV).
 * Some modifications have been made. The original note is reproduced here:
 *
 * This source file is used to print out a stack-trace when your program
 * segfaults. It is relatively reliable and spot-on accurate.
 *
 * This code is in the public domain. Use it as you see fit, some credit
 * would be appreciated, but is not a prerequisite for usage. Feedback
 * on it's use would encourage further development and maintenance.
 *
 * Due to a bug in gcc-4.x.x you currently have to compile as C++ if you want
 * demangling to work.
 *
 * Please note that it's been ported into my ULS library, thus the check for
 * HAS_ULSLIB and the use of the sigsegv_outp macro based on that define.
 *
 * Author: Jaco Kroon <jaco@kroon.co.za>
 *
 * Copyright (C) 2005 - 2010 Jaco Kroon
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/* Bug in gcc prevents from using CPP_DEMANGLE in pure "C" */
#if !defined(__cplusplus) && !defined(NO_CPP_DEMANGLE)
#define NO_CPP_DEMANGLE
#endif

#include <openlat/config.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#ifndef NO_CPP_DEMANGLE
#include <cxxabi.h>
#ifdef __cplusplus
using __cxxabiv1::__cxa_demangle;
#endif
#endif



#ifdef HAVE_EXECINFO_H

#include <ucontext.h>
#include <execinfo.h>
#include <dlfcn.h>

void dump_backtrace() {
  Dl_info mdlinfo;
  char syscom[256];
  int f = 0;
  char funcname[1024];
  char fileline[1024];
  int status;

  void *_bt[100];
  void **bt = (void **)_bt;
  int sz = backtrace(bt, 100);
  // skip i = 0 since it is this dump_backtrace function
  for(int i = 1; i < sz; ++i) {
    if(!dladdr(bt[i], &mdlinfo)) break;
    if (mdlinfo.dli_saddr == NULL) continue;

    const char *symname = mdlinfo.dli_sname;

#ifndef NO_CPP_DEMANGLE
    char * tmp = __cxa_demangle(symname, NULL, 0, &status);

    if (status == 0 && tmp) symname = tmp;
#endif

#ifdef ADDR2LINEBIN
    sprintf(syscom, ADDR2LINEBIN " --demangle --basenames --functions -e %s %p", mdlinfo.dli_fname, mdlinfo.dli_saddr);
    FILE *cmd = popen(syscom, "r");
    int num = fscanf(cmd, "%s\n%s\n", funcname, fileline);
    status = pclose(cmd);
    if (num == EOF || WIFEXITED(status)) {
      status = WEXITSTATUS(status);
    } 
    if (status == 0 && strcmp(funcname, "??") != 0) {
      symname = funcname;
    }
    else {
      sprintf(fileline, "in %s", mdlinfo.dli_fname);
    }
#else
    sprintf(fileline, "in %s", mdlinfo.dli_fname);
#endif

    if (++f == 1) {
      fprintf(stderr, "   at");
    }
    else {
      fprintf(stderr, "   by");
    }

    fprintf(stderr, " %p: %s (%s)\n", mdlinfo.dli_saddr, symname, fileline);


#ifndef NO_CPP_DEMANGLE
    if (tmp) free(tmp);
#endif

    if(mdlinfo.dli_sname && !strcmp(mdlinfo.dli_sname, "main")) break;

  }
}

static void signal_segv(int signum, siginfo_t* info, void* ptr) {
  signal(SIGSEGV, SIG_DFL);

  //static const char *si_codes[3] = {"", "SEGV_MAPERR", "SEGV_ACCERR"};
  static const char *si_message[3] = {"", "Address not mapped", "Invalid permissions"};

  fprintf(stderr, "Segmentation Fault!\n\n");
  fprintf(stderr, "%s\n", si_message[info->si_code]);

  dump_backtrace();

  char msg[256]; 
  sprintf(msg, " Address %p cannot be accessed", info->si_addr);
  if (info->si_errno != 0) perror(msg);
  else fprintf(stderr, "%s\n", msg);
  //_exit (-1);
  // avoid unused warning
  signum = 0; ptr = NULL;
}

//XXX: we have removed the constructor attribute so that initialization
// has to be made explicitely
// void __attribute__((constructor)) setup_sigsegv() {
void setup_sigsegv() {
  struct sigaction action;
  memset(&action, 0, sizeof(action));
  action.sa_sigaction = signal_segv;
  action.sa_flags = SA_SIGINFO;
  if(sigaction(SIGSEGV, &action, NULL) < 0)
    perror("sigaction");
}

#else
void dump_backtrace() {
}

void setup_sigsegv() {
}
#endif //HAVE_BACKTRACE
