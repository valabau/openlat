%pure-parser
%name-prefix="htk_"
%locations
%defines
%debug
%error-verbose
%parse-param { openlat::HtkContext* context }
%lex-param { void* scanner  }

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
%token <ival> XSCALE XSCORE
%token <dval> FLOAT
%token <sval> STRING
%token <sval> UNK_OPTION
%token        SLF_VERSION UTTERANCE SUBLAT BASE LMNAME LMSCALE WDPENALTY NDPENALTY ECSCALE
%token        LMINSCALE LMOUTSCALE WDPENALTY_OUTPUT LMIN LMOUT OUTPUT
%token        ACSCALE AMSCALE VOCAB COVERAGE
%token        INITIAL_NODE FINAL_NODE FEATURES
%token        NODES LINKS
%token        NODE TIME WORD SUBS VAR
%token        LINK START_NODE END_NODE DIV ACOUSTIC NGRAM LANGUAGE POSTERIOR
%token        EXPECTED_ERROR CORRECT ESTIMATED_CORRECT
%token        ERR

/* Type definition for rules */
%type  <dval> number

/* %printer    { debug_stream () << *$$; } STRING UNK_OPTION */
%destructor { delete[] $$; } STRING UNK_OPTION

/* %printer    { debug_stream () << $$; } INT FLOAT */


%{
	#include <iostream>
	#include "htk.h"
	#include <openlat/utils.h>

	using namespace std;
	using namespace openlat;

	int htk_lex(YYSTYPE* lvalp, YYLTYPE* llocp, void* scanner);		

	void htk_error(YYLTYPE* locp, openlat::HtkContext* context, const char* err)
	{
		cout << locp->first_line << ":" << err << endl;
	}

	#define scanner context->scanner
	
%}

%%

/*

The following rules define the syntax  of a lattice using the first character of
the field names defined below to represent the entire field. Any unrecognised
fields may be ignored and no user defined fields may share first characters with
pre-defined field names. parentheses () indicate that terms within may appear
in any order and braces {} indicate that terms within are optional.

lattice_file = header size_def terms
header       = ( V {\n} U {\n} S {\n} ) { base {\n} start {\n} end {\n} dir {\n}
                                   tscale {\n} lmscale {\n} lmpen {\n} }
size_def     = ( L N ) \n
terms        = node_def | link_def { terms }
node_def     = I { t W v d a } \n
link_def     = J ( S E { W v d a n l} ) \n

*/

%start lattice_file;

/* File structure */
lattice_file: { /* Init variables */ }
              header
              { /* Check compulsory paramenters */
              
                // Make this first to ensure that the initial state is state 0
              	if (not context->htk.start_name.empty()) {
              	  context->htk.getStateId(context->htk.start_name);
              	} 
              }
              size_def
              { /* Check valid lattice size */ }
              terms
              { }

/* header options */
header: /* empty */
  | header option

option: SLF_VERSION STRING { context->htk.version = string($2); delete[] $2; }   ENDL
      | UTTERANCE   STRING { context->htk.utterance = string($2); delete[] $2; } ENDL
      | SUBLAT      STRING { context->htk.sublat = string($2); delete[] $2; }    ENDL
      | BASE        number { context->htk.base = $2; }
      | LMNAME      STRING { context->htk.lmname = string($2); delete[] $2; }    ENDL
      | LMSCALE     number { context->htk.assignWeight(HtkLattice::LANGUAGE, $2); }
      | LMINSCALE   number { context->htk.assignWeight(HtkLattice::LMIN, $2); }
      | LMOUTSCALE  number { context->htk.assignWeight(HtkLattice::LMOUT, $2); }
      /* Word Penalty in log, as SRILM does */
      //| WDPENALTY   number { context->htk.weights.wdpenalty = context->htk.changeBase($2); }
      | WDPENALTY   number { context->htk.assignWeight(HtkLattice::WDPENALTY, $2); }
      | NDPENALTY   number { context->htk.assignWeight(HtkLattice::NOISE_PENALTY, $2); }
      | WDPENALTY_OUTPUT number { context->htk.assignWeight(HtkLattice::OUTPUT_WDPENALTY, $2); }
      | ACSCALE     number { context->htk.assignWeight(HtkLattice::ACOUSTIC, $2); }
      | AMSCALE     number { context->htk.weights.amscale = $2; }
      | VOCAB       STRING { context->htk.vocab = string($2); delete[] $2; }   ENDL
      | XSCALE      number { context->htk.assignWeight(HtkLattice::feature_t(HtkLattice::XSCORE + $1-1), $2); }
      | INITIAL_NODE  INT  { context->htk.start_name = to_string($2); }
      | FINAL_NODE    INT  { context->htk.end_name   = to_string($2); }
      | UNK_OPTION  STRING ENDL { delete[] $1; delete[] $2; }
      | ENDL

/* Size definitions */
size_def: /* empty */
        | size_def def

def: NODES INT { context->num_nodes = $2; }
   | LINKS INT { context->num_links = $2; }
   | ENDL


/* Nodes and links */
terms: /* empty */
     | terms term

term: node
    | link
    | ENDL

/* Node */
node: NODE INT
      { 
        context->node_ptr = context->htk.getState(to_string($2));
        context->link_ptr = &context->node_ptr->link; 
      }
      node_options
      { /* Post-process node */
        context->node_ptr = 0;
        context->link_ptr = 0;     
      }
      ENDL
node_options: /* empty */
            | node_options node_option
node_option: TIME number { context->node_ptr->time = $2; }
           | SUBS STRING { context->node_ptr->sublat = string($2); delete[] $2; }
           | COVERAGE  STRING { delete[] $2; }
           | link_option

           /* Link info in node   
           | VAR  INT    { context->link_ptr->var = $2; }
           | ACOUSTIC   number { context->htk.assign(context->link_ptr, HtkLattice::ACOUSTIC, $2); }
           | DIV        STRING { context->link_ptr->div = string($2); delete[] $2; }
           | WORD       STRING { context->htk.assignInput(context->link_ptr, string($2)); delete[] $2; }
           | OUTPUT     STRING { context->htk.assignOutput(context->link_ptr, string($2)); delete[] $2; }
           | UNK_OPTION  STRING { delete[] $1; delete[] $2; }
           */

/* Links */
link: LINK INT
      {
        context->link_ptr = context->htk.getArc(to_string($2));
      }
      link_options
      {
        context->link_ptr = NULL;
      }
      ENDL
link_options: /* empty */
            | link_options link_option
link_option: START_NODE INT    { context->link_ptr->start = context->htk.getStateId(to_string($2)); }
           | END_NODE   INT    { context->link_ptr->end   = context->htk.getStateId(to_string($2)); }
           | WORD       STRING { context->htk.assignInput(context->link_ptr, string($2)); delete[] $2; }
           | OUTPUT     STRING { context->htk.assignOutput(context->link_ptr, string($2)); delete[] $2; }
           | VAR        INT    { context->link_ptr->var = $2; }
           | DIV        STRING { context->link_ptr->div = string($2); delete[] $2; }
           | ACOUSTIC   number { context->htk.assign(context->link_ptr, HtkLattice::ACOUSTIC, $2); }
           | LANGUAGE   number { context->htk.assign(context->link_ptr, HtkLattice::LANGUAGE, $2); }
           | NGRAM      number { context->htk.assign(context->link_ptr, HtkLattice::NGRAM, $2); }
           | LMIN       number { context->htk.assign(context->link_ptr, HtkLattice::LMIN, $2); }
           | LMOUT      number { context->htk.assign(context->link_ptr, HtkLattice::LMOUT, $2); }
           | POSTERIOR  number { context->htk.assign(context->link_ptr, HtkLattice::POSTERIOR, $2); }
           | XSCORE     number { context->htk.assign(context->link_ptr, HtkLattice::feature_t(HtkLattice::XSCORE + $1-1), $2); }
           | UNK_OPTION STRING { delete[] $1; delete[] $2; }

/* Helper rules */
number:	FLOAT { $$ = $1; }
      | INT   { $$ = (double) $1; }

