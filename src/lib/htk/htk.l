/* based on http://www.phpcompiler.org/articles/reentrantparser.html */
%option reentrant
%option prefix="openlat_htk_"
%option bison-bridge
%option bison-locations
%option noyywrap
%option yylineno
%option batch 
%option debug 
%option case-insensitive

%{

#include "htk.h"

using namespace openlat;
using namespace openlat::htk;


/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return parser::token::END


  #define MAX_COND_LEVEL 128
  #undef  BEGIN
  #define BEGIN yyextra->stack.push_back(yyg->yy_start); (yyg->yy_start) = 1 + 2 *
  #define END() (yyg->yy_start) = yyextra->stack.back(); yyextra->stack.pop_back(); 


%}

%x INCOMMENT LONGSTR NODELINE LINKLINE NWSTR

WSCHAR              [ \t\f\v\r]
NOT_WSCHAR          [^ \t\f\v\r\n]
OCTDIGIT            [0-7]
HEXDIGIT            [0-9a-fA-F]
DIGIT               [0-9]

                /* Strings */
ESCAPE              (\\[ntvbrfa\n\\?'"])
IN_STR              ([^\\\n]|{ESCAPE})*
/* IN_STR              ([^\\\n"]|{ESCAPE})* */
STRING              \"{IN_STR}\"
CHARLIT             \'{IN_STR}\'

                /* Intergers */
OCTAL               (0{OCTDIGIT}+)
HEXADECIMAL         (0[xX]{HEXDIGIT}+)
UINTEGER            ([1-9]{DIGIT}*)|0
INTEGER             [+-]?{UINTEGER}

                /* Floats */
EXPONENT            [Ee][+-]?{DIGIT}+
FRACTIONAL          ({DIGIT}+\.)|({DIGIT}*\.{DIGIT}+)
FLOAT               [+-]?{FRACTIONAL}{EXPONENT}?|[+-]?{DIGIT}+{EXPONENT}|-inf

                /* Others */
WHITESPACE          {WSCHAR}+|({WSCHAR}*\\)+{WSCHAR}*
NEWLINE             [\n]

IDENTCHARS	    [_/a-zA-Z]
IDENT               {IDENTCHARS}({IDENTCHARS}|[0-9])*


%%

"#".*                   { /* yylloc->lines (1); ylloc->step (); return(COMMENT);*/ }

      /* Header */
("VERSION"|"V"){WHITESPACE}*"="   { BEGIN(LONGSTR); return parser::token::SLF_VERSION; }
("UTTERANCE"|"U"){WHITESPACE}*"=" { BEGIN(LONGSTR); return parser::token::UTTERANCE; }
("SUBLAT"|"S"){WHITESPACE}*"="    { BEGIN(LONGSTR); return parser::token::SUBLAT; }
"base"{WHITESPACE}*"="            { return parser::token::BASE; }
"lmname"{WHITESPACE}*"="          { BEGIN(LONGSTR); return parser::token::LMNAME; }
"lmscale"{WHITESPACE}*"="         { return parser::token::LMSCALE; }
"lminscale"{WHITESPACE}*"="       { return parser::token::LMINSCALE; }
"lmoutscale"{WHITESPACE}*"="      { return parser::token::LMOUTSCALE; }
"wdpenalty_output"{WHITESPACE}*"=" { return parser::token::WDPENALTY_OUTPUT; }
"wdpenalty"{WHITESPACE}*"="       { return parser::token::WDPENALTY; }
"ndpenalty"{WHITESPACE}*"="       { return parser::token::NDPENALTY; }
"acscale"{WHITESPACE}*"="         { return parser::token::ACSCALE; }
"amscale"{WHITESPACE}*"="         { return parser::token::AMSCALE; }
"start"{WHITESPACE}*"="           { return parser::token::INITIAL_NODE; }
"end"{WHITESPACE}*"="             { return parser::token::FINAL_NODE; }
"x"{INTEGER}"scale"{WHITESPACE}*"=" { errno = 0;
                                    long n = strtol(yytext + 1,NULL,10); 
                                    if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                                      htk_error(yylloc, yyextra, "integer is out of range");
                                    }
                                    yylval->ival = n;
                         
                                    return parser::token::XSCALE; 
                                  }
"vocab"{WHITESPACE}*"="           { BEGIN(LONGSTR); return parser::token::VOCAB; }

<LONGSTR>{
  {IN_STR}     { // Copy and don't remove "" and '' 
                 yylval->sval = parser_txtcpy(yytext, false);
                 return parser::token::STRING; 
               }
  {NEWLINE}    { 
                 END(); 
                 // yylloc->lines (yyleng); yylloc->step (); 
                 return parser::token::ENDL;
               }
}


      /* Size defs */
("NODES"|"N"){WHITESPACE}*"="   { return parser::token::NODES; }
("LINKS"|"L"){WHITESPACE}*"="   { return parser::token::LINKS; }

      /* Nodes */
"I"{WHITESPACE}*"="             { BEGIN(NODELINE); return parser::token::NODE; }
<NODELINE>{
  ("time"|"t"){WHITESPACE}*"="   { return parser::token::TIME; }
  ("WORD"|"W"){WHITESPACE}*"="   { BEGIN(NWSTR); return parser::token::WORD; }
  "L"{WHITESPACE}*"="            { return parser::token::SUBS; }
  ("var"|"v"){WHITESPACE}*"="    { return parser::token::VAR; }
  ("div"|"d"){WHITESPACE}*"="      { return parser::token::DIV; }
  ("acoustic"|"a"){WHITESPACE}*"=" { return parser::token::ACOUSTIC; }
  ("coverage"|"c"){WHITESPACE}*"=" { BEGIN(NWSTR); return parser::token::COVERAGE; }
  ("output"|"O"){WHITESPACE}*"="   { BEGIN(NWSTR); return parser::token::OUTPUT; }
  "#".*                            { /* yylloc->lines (1); ylloc->step (); return(parser::token::COMMENT);*/ }
  {WHITESPACE}                   { /* yylloc->step(); */ }
  {NEWLINE}                      { 
                                   END(); 
                                   /* yylloc->lines (yyleng); yylloc->step (); */
                                   return parser::token::ENDL;
                                 }
}

      /* Links */
"J"{WHITESPACE}*"="                { BEGIN(LINKLINE); return parser::token::LINK; }
<LINKLINE>{
  ("START"|"S"){WHITESPACE}*"="    { return parser::token::START_NODE; }
  ("END"|"E"){WHITESPACE}*"="      { return parser::token::END_NODE; }
  ("WORD"|"W"){WHITESPACE}*"="     { BEGIN(NWSTR); return parser::token::WORD; }
  ("div"|"d"){WHITESPACE}*"="      { return parser::token::DIV; }
  ("acoustic"|"a"){WHITESPACE}*"=" { return parser::token::ACOUSTIC; }
  ("ngram"|"n"){WHITESPACE}*"="    { return parser::token::NGRAM; }
  ("language"|"l"){WHITESPACE}*"=" { return parser::token::LANGUAGE; }
  ("lmin"){WHITESPACE}*"="         { return parser::token::LMIN; }
  ("lmout"){WHITESPACE}*"="        { return parser::token::LMOUT; }
  ("posterior"|"p"){WHITESPACE}*"=" { return parser::token::POSTERIOR; }
  "x"{INTEGER}{WHITESPACE}*"="     { errno = 0;
                                    long n = strtol(yytext + 1,NULL,10); 
                                    if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                                      htk_error(yylloc, yyextra, "integer is out of range");
                                    }
                                    yylval->ival = n;
                         
                                    return parser::token::XSCORE; 
                                  }
  ("output"|"O"){WHITESPACE}*"="    { BEGIN(NWSTR); return parser::token::OUTPUT; }
  "#".*                            { /* yylloc->lines (1); ylloc->step (); return(Cparser::token::OMMENT);*/ }
  {WHITESPACE}                     { /* yylloc->step(); */ }
  {NEWLINE}                        { END(); 
                                     // yylloc->lines (yyleng); yylloc->step (); 
                                     return parser::token::ENDL;
                                   }
}


 
{NEWLINE}       { 
                  // yylloc->lines (yyleng); yylloc->step (); 
                  return parser::token::ENDL;
                }

{WHITESPACE}    { /* yylloc->step(); *//*return(parser::token::SPACE);*/ }

<NWSTR>{
  {NEWLINE}     { 
                  yylval->sval = parser_txtcpy("", false);
                  unput('\n'); 
                  END();
                  return parser::token::STRING; 
                }

  {WHITESPACE}  { 
                  yylval->sval = parser_txtcpy("", false); 
                  END();
                  return parser::token::STRING; 
                }

  {NOT_WSCHAR}* {
                  yylval->sval = parser_txtcpy(yytext, false); 
                  END();
                  return parser::token::STRING; 
                }
}


<*>{

  {OCTAL}	{ 
  		        errno = 0;
  	          long n = strtol(yytext+1,NULL,8); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                htk_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
              return parser::token::INT;
  		}
  
  {HEXADECIMAL}   {
  	          errno = 0;
  	          long n = strtol(yytext+2,NULL,16); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                htk_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
  	          return parser::token::INT;
  		}
  
  {INTEGER}       { 
  	          errno = 0;
  	          long n = strtol(yytext,NULL,10); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                htk_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
              return parser::token::INT;
  		}
  
  
  {FLOAT}         { 
  	          errno = 0;
  	          double n = strtod(yytext,NULL); 
              if (! (errno != ERANGE)) {
                htk_error(yylloc, yyextra, "float is out of range");
              }
  	          yylval->dval = n;
              return parser::token::FLOAT;
  		}
}

<INITIAL,NODELINE,LINKLINE>{STRING} |
{CHARLIT}       { // Copy and remove "" and '' 
                  yylval->sval = parser_txtcpy(yytext, true); 
                  return parser::token::STRING; 
                }


<*>{IDENT}{WHITESPACE}*"=" { 
                  yylval->sval = parser_txtcpy(yytext, false); 
                  BEGIN(NWSTR); 
                  return parser::token::UNK_OPTION; 
                }

<*>.            { htk_error(yylloc, yyextra, "invalid character"); exit(-1); }

%%

void openlat::htk::HtkContext::init_scanner()
{
	yylex_init(&scanner);
	yyset_extra(this, scanner);
}

void openlat::htk::HtkContext::destroy_scanner()
{
	yylex_destroy(scanner);
  	scanner = 0;
}

void openlat::htk::htk_error(YYLTYPE* locp, openlat::htk::HtkContext*, const char* err) 
{
	openlat::parser_error<YYLTYPE>(*locp, string(err));
}

void openlat::htk::parser::error(const YYLTYPE &locp, const string &err) 
{
	openlat::parser_error<YYLTYPE>(locp, err);
}