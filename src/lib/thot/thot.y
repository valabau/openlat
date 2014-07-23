// based on http://www.phpcompiler.org/articles/reentrantparser.html
%pure-parser
%skeleton "lalr1.cc"                          /*  -*- C++ -*- */
%define namespace "openlat::thot"
%name-prefix="openlat_thot_"
%locations
%defines
%debug
%error-verbose
%parse-param { openlat::thot::ThotContext* context }
%lex-param { void* scanner  }

%{
	#include <iostream>
	#include "thot.h"
%}

// Symbols.
%union
{
  int          ival;
  double       dval;
  char	      *sval;
};

/* Tokens from scanner */
%token        END 0 "end of file"
%token        ENDL  "end of line"
%token <ival> INT
%token <dval> FLOAT
%token <sval> STRING

/* Type definition for rules */
%type  <dval> number
%type  <sval> words

/* %printer    { debug_stream () << *$$; } STRING */
%destructor { delete[] $$; } STRING words

/* %printer    { debug_stream () << $$; } INT FLOAT */


%{
	#include <openlat/utils.h>

	using namespace std;
	using namespace openlat;

	int openlat_thot_lex(YYSTYPE* lvalp, YYLTYPE* llocp, void* scanner);		

	#define scanner context->scanner
	
%}

%%

/*

The following rules define the syntax  of a thot word graph 

wgp_file = states arcs
states   = id+ \n
arcs     = start_node end_node score words \n
words    = word+
id, start_node, end_node <- int
word <- string

*/

%start wgp_file;

/* File structure */
wgp_file:     { /* Init variables */ }
			  empty_lines
              states ENDL 
			  empty_lines
              { cerr << "Start reading arcs\n"; }
              arcs
              { /* Check consistency */}

/* header options */
empty_lines: /* empty */ | ENDL | empty_lines ENDL
states: /* empty */
  | states state

/* State */
state: INT
      { 
        context->thot.getState(to_string($1));
      }
  

/* Nodes and links */
arcs: /* empty */
     | arcs arc

/* Links */
arc: INT INT number words
      {
        context->link_ptr = context->thot.getNewArc();
        context->link_ptr->start = context->thot.getStateId(to_string($1));
        context->link_ptr->end   = context->thot.getStateId(to_string($2));
        context->link_ptr->score = $3;
        context->thot.assignWords(context->link_ptr, string($4)); delete[] $4;
        context->link_ptr = NULL;
      }
      ENDL
    | ENDL
      
/* Helper rules */

words: { $$ = 0; } /* empty */
	 | words STRING {
	   $$ = new char[strlen($1)+strlen($2)+1];
	   strcpy($$, $1);
	   delete[] $1; 	
	   strcat($$, $2);
	   delete[] $2; 	
     }

number:	FLOAT { $$ = $1; }
      | INT   { $$ = (double) $1; }

