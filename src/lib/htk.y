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
	#include <openlat/htk-compiler.h>

	using namespace std;

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
lattice_file: { /* Init variables */
                node_ptr = NULL;
                link_ptr = NULL;
                if (logbase != FLOAT_NONE) driver.slf.base = logbase;
                if (driver.slf.features.empty()) load_features = true;
              }
              header
              { /* Check compulsory paramenters */
                /* Fix feature mismatches*/
                if (driver.slf.scorer.feature_scales.size() > 0) {
                  if (driver.slf.features.empty()) {
                    size_t f = 0;
                    stringstream ss;
                    ss << "f" << f;
                    for (f = 1; f < driver.slf.scorer.feature_scales.size(); f++) {
                      ss << ",f" << f;
                    }
                    driver.slf.setFeatures(ss.str());
                  }
                  else if (driver.slf.scorer.feature_scales.size() < driver.slf.features.size()) {
                    driver.slf.scorer.feature_scales.resize(driver.slf.features.size(), 1.0);
                  }
                }
                assert_bt(driver.slf.features.size() == driver.slf.scorer.feature_scales.size(), "Mismatch in the number of features and feature scores");
              }
              size_def
              { /* Check valid lattice size */}
              terms
              {
                if (driver.slf.initial != UNK_NODE) {
                  driver.slf.initial = driver.slf.getInternalNodeId(driver.slf.initial);
                }
                if (driver.slf.final != UNK_NODE) {
                  driver.slf.final = driver.slf.getInternalNodeId(driver.slf.final);
                }
              }

/* header options */
header: /* empty */
  | header option

option: SLF_VERSION STRING { driver.slf.version = string($2); delete[] $2; }   ENDL
      | UTTERANCE   STRING { driver.slf.utterance = string($2); delete[] $2; } ENDL
      | SUBLAT      STRING { driver.slf.sublat = string($2); delete[] $2; }    ENDL
      | BASE        number { if (logbase == FLOAT_NONE) driver.slf.base = $2; }
      | LMNAME      STRING { driver.slf.scorer.lmname = string($2); delete[] $2; }    ENDL
      | LMSCALE     number { driver.slf.scorer.lmscale = $2; }
      | LMINSCALE   number { driver.slf.scorer.lminscale = $2; }
      | LMOUTSCALE  number { driver.slf.scorer.lmoutscale = $2; }
      /* Word Penalty in log, as SRILM does */
      //| WDPENALTY   number { driver.slf.scorer.wdpenalty = driver.slf.changeBase($2); }
      | WDPENALTY   number { driver.slf.scorer.wdpenalty = $2; }
      | NDPENALTY   number { driver.slf.scorer.ndpenalty = $2; }
      //| WDPENALTY_OUTPUT number { driver.slf.scorer.wdpenalty_output = driver.slf.changeBase($2); }
      | WDPENALTY_OUTPUT number { driver.slf.scorer.wdpenalty_output = $2; }
      | ACSCALE     number { driver.slf.scorer.acscale = $2; }
      | AMSCALE     number { driver.slf.scorer.amscale = $2; }
      | ECSCALE     number { driver.slf.scorer.ecscale = $2; }
      | VOCAB       STRING { driver.slf.vocab = $2; delete[] $2; }   ENDL
      | FEATURES    STRING { if (load_features) driver.slf.setFeatures($2); delete[] $2; }   ENDL
      | XSCALE      number { if (load_features) driver.slf.setFeatureScale($1-1, $2); }
      | INITIAL_NODE  INT  { driver.slf.initial = static_cast<node_id_t>($2); }
      | FINAL_NODE    INT  { driver.slf.final = static_cast<node_id_t>($2); }
      | UNK_OPTION  STRING ENDL { delete[] $1; delete[] $2; }
      | ENDL

/* Size definitions */
size_def: /* empty */
        | size_def def

def: NODES INT { num_nodes = $2; }
   | LINKS INT { num_links = $2; }
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
        node_ptr = new Node();
        node_ptr->id = static_cast<node_id_t>($2);
        node_ptr->original_id = static_cast<node_id_t>($2);
      }
      node_options
      {
        driver.slf.addNode(node_ptr);
        node_ptr = NULL;
      }
      ENDL
node_options: /* empty */
            | node_options node_option
node_option: TIME number { node_ptr->time = $2; }
           | SUBS STRING { node_ptr->sublat = string($2); delete[] $2; }
           | VAR  INT    { node_ptr->var = $2; }
           | COVERAGE  STRING { delete[] $2; }   
           | ACOUSTIC   number { node_ptr->acoustic = driver.slf.changeBase($2); }
           | DIV        STRING { node_ptr->div = string($2); delete[] $2; }
           | WORD       STRING { 
                                 vector<string> inwords;
                                 tokenize($2, inwords, "_");
                                 for (size_t w = 0; w < inwords.size(); w++) {
                                   node_ptr->input.push_back(Word(inwords[w]));
                                 }
                                 delete[] $2;
                               }
           | OUTPUT     STRING {
                                 vector<string> outwords;
                                 tokenize($2, outwords, "_");
                                 for (size_t w = 0; w < outwords.size(); w++) {
                                   node_ptr->output.push_back(Word(outwords[w]));
                                 }
                                 delete[] $2;
                               }
           | UNK_OPTION  STRING { delete[] $1; delete[] $2; }

/* Links */
link: LINK INT
      {
        link_ptr = new Link(driver.slf.features.size());
        fill(link_ptr->features, link_ptr->features + driver.slf.features.size(), .0);
        link_ptr->original_id = static_cast<link_id_t>($2);
        has_word = false;
        has_ac = false;
      }
      link_options
      {
        if (link_ptr->start != UNK_NODE and link_ptr->end != UNK_NODE) {
          std::ostringstream ss;
          if (not has_word) {
            link_ptr->input = driver.slf.nodes[link_ptr->end]->input;
            link_ptr->output = driver.slf.nodes[link_ptr->end]->output;
          }
          if (not has_ac) link_ptr->acoustic = driver.slf.nodes[link_ptr->end]->acoustic; 
          
          if (link_ptr->input.empty()) {
            ss << WordFactory::null;
          }
          else {
            ss << link_ptr->input[0];
            for (size_t i = 1; i < link_ptr->input.size(); i++) {
              ss << "_" << link_ptr->input[i];
            }
          }
          if (driver.slf.has_joint_phrases and not link_ptr->output.empty()) {
            ss << "|||" << link_ptr->output[0];
            for (size_t o = 1; o < link_ptr->output.size(); o++) {
              ss << "_" << link_ptr->output[o];
            }
          }
          link_ptr->word_id = driver.slf.wordmap.getWordId(Word(ss.str()));

          driver.slf.addLink(link_ptr);
        }
        else {
          cerr << "Warning: Discarding link with non-existing start or end nodes\n";
          delete link_ptr;
        }
        link_ptr = NULL;
      }
      ENDL
link_options: /* empty */
            | link_options link_option
link_option: START_NODE INT    { link_ptr->start = driver.slf.getInternalNodeId($2); }
           | END_NODE   INT    { link_ptr->end = driver.slf.getInternalNodeId($2); }
           | WORD       STRING {
                                 has_word = true;
                                 vector<string> inwords;
                                 tokenize($2, inwords, "_");
                                 for (size_t w = 0; w < inwords.size(); w++) {
                                   link_ptr->input.push_back(Word(inwords[w]));
                                 }
                                 delete[] $2;
                               }
           | VAR        INT    { link_ptr->var = $2; }
           | DIV        STRING { link_ptr->div = string($2); delete[] $2; }
           | ACOUSTIC   number { link_ptr->acoustic = driver.slf.changeBase($2); has_ac = true; }
           | NGRAM      number { link_ptr->ngram = driver.slf.changeBase($2); }
           | LANGUAGE   number { link_ptr->language = driver.slf.changeBase($2); }
           | LMIN       number { link_ptr->language_in = driver.slf.changeBase($2); }
           | LMOUT      number { link_ptr->language_out = driver.slf.changeBase($2); }
           | POSTERIOR  number { link_ptr->posterior = driver.slf.changeBase($2); }
           | EXPECTED_ERROR number { link_ptr->error = driver.slf.changeBase($2); }
           | CORRECT    number { link_ptr->correct = bool($2); driver.slf.tagged = true; }
           | ESTIMATED_CORRECT    number { link_ptr->estimated_correct = bool($2); }
           | FEATURES   STRING { if (load_features) link_ptr->setFeatures($2); delete[] $2; }
           | XSCORE     number { if (load_features) {
                                   assert_bt($1-1 < static_cast<int>(driver.slf.features.size()), 
                                     "Invalid feature index "<<($1-1)<<" of "<<driver.slf.features.size()<<" features");
                                   link_ptr->features[$1-1] = $2; 
                                 }
                               }
           | OUTPUT     STRING {
                                 has_word = true;
                                 vector<string> outwords;
                                 tokenize($2, outwords, "_");
                                 for (size_t w = 0; w < outwords.size(); w++) {
                                   link_ptr->output.push_back(Word(outwords[w]));
                                 }
                                 delete[] $2;
                               }
           | UNK_OPTION STRING { delete[] $1; delete[] $2; }

/* Helper rules */
number:	FLOAT { $$ = $1; }
      | INT   { $$ = (double) $1; }

