%option reentrant
%option prefix="htk_"
%option bison-bridge
%option bison-locations
%option noyywrap
%option yylineno
%option batch 
%option debug 
%option case-insensitive

%{
	#include "htk.h"
	#include "yacc.htk.hpp"

  using namespace openlat;

  /* Utilities */
  char *htk_txtcpy(const char *txt, bool remove_quotes) {
    char *ptr;
    int len = (remove_quotes)? strlen(txt)-2:strlen(txt);
    ptr = new char[len+1];
    strncpy(ptr,txt + ((remove_quotes)?1:0), len);
    ptr[len] = '\0';
    return ptr;
  }
  
	void htk_error(YYLTYPE* locp, HtkContext* context, const char* err);

  #define MAX_COND_LEVEL 128
  #undef  BEGIN
  #define BEGIN htk_cond_stack[htk_cond_idx++] = yyg->yy_start; (yyg->yy_start) = 1 + 2 *
  #define END() (yyg->yy_start) = htk_cond_stack[--htk_cond_idx] 
  int htk_cond_idx = 0;
  int htk_cond_stack[MAX_COND_LEVEL];


	#define YY_EXTRA_TYPE HtkContext*
	#define YY_USER_ACTION yylloc->first_line = yylineno;

	#define YY_INPUT(buf,result,max_size) \
	{ \
		char c; \
		yyextra->is >> c; \
		if (yyextra->is.eof()) \
			result = YY_NULL; \
		else { \
			buf[0] = c; \
			result = 1; \
		} \
	}
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
("VERSION"|"V"){WHITESPACE}*"="   { BEGIN(LONGSTR); return SLF_VERSION; }
("UTTERANCE"|"U"){WHITESPACE}*"=" { BEGIN(LONGSTR); return UTTERANCE; }
("SUBLAT"|"S"){WHITESPACE}*"="    { BEGIN(LONGSTR); return SUBLAT; }
"base"{WHITESPACE}*"="            { return BASE; }
"lmname"{WHITESPACE}*"="          { BEGIN(LONGSTR); return LMNAME; }
"lmscale"{WHITESPACE}*"="         { return LMSCALE; }
"lminscale"{WHITESPACE}*"="       { return LMINSCALE; }
"lmoutscale"{WHITESPACE}*"="      { return LMOUTSCALE; }
"wdpenalty_output"{WHITESPACE}*"=" { return WDPENALTY_OUTPUT; }
"wdpenalty"{WHITESPACE}*"="       { return WDPENALTY; }
"ndpenalty"{WHITESPACE}*"="       { return NDPENALTY; }
"acscale"{WHITESPACE}*"="         { return ACSCALE; }
"amscale"{WHITESPACE}*"="         { return AMSCALE; }
"ecscale"{WHITESPACE}*"="         { return ECSCALE; }
"start"{WHITESPACE}*"="           { return INITIAL_NODE; }
"end"{WHITESPACE}*"="             { return FINAL_NODE; }
"x"{INTEGER}"scale"{WHITESPACE}*"=" { errno = 0;
                                    long n = strtol(yytext + 1,NULL,10); 
                                    if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                                      htk_error(yylloc, yyextra, "integer is out of range");
                                    }
                                    yylval->ival = n;
                         
                                    return XSCALE; 
                                  }
("features"|"f"){WHITESPACE}*"="  { BEGIN(NWSTR); return FEATURES; }
"vocab"{WHITESPACE}*"="           { BEGIN(LONGSTR); return VOCAB; }

<LONGSTR>{
  {IN_STR}     { // Copy and don't remove "" and '' 
                 yylval->sval = htk_txtcpy(yytext, false); 
                 return STRING; 
               }
  {NEWLINE}    { 
                 END(); 
                 // yylloc->lines (yyleng); yylloc->step (); 
                 return ENDL;
               }
}


      /* Size defs */
("NODES"|"N"){WHITESPACE}*"="   { return NODES; }
("LINKS"|"L"){WHITESPACE}*"="   { return LINKS; }

      /* Nodes */
"I"{WHITESPACE}*"="             { BEGIN(NODELINE); return NODE; }
<NODELINE>{
  ("time"|"t"){WHITESPACE}*"="   { return TIME; }
  ("WORD"|"W"){WHITESPACE}*"="   { BEGIN(NWSTR); return WORD; }
  "L"{WHITESPACE}*"="            { return SUBS; }
  ("var"|"v"){WHITESPACE}*"="    { return VAR; }
  ("div"|"d"){WHITESPACE}*"="      { return DIV; }
  ("acoustic"|"a"){WHITESPACE}*"=" { return ACOUSTIC; }
  ("coverage"|"c"){WHITESPACE}*"=" { BEGIN(NWSTR); return COVERAGE; }
  ("output"|"O"){WHITESPACE}*"="   { BEGIN(NWSTR); return OUTPUT; }
  "#".*                            { /* yylloc->lines (1); ylloc->step (); return(COMMENT);*/ }
  {WHITESPACE}                   { /* yylloc->step(); */ }
  {NEWLINE}                      { 
                                   END(); 
                                   /* yylloc->lines (yyleng); yylloc->step (); */
                                   return ENDL;
                                 }
}

      /* Links */
"J"{WHITESPACE}*"="                { BEGIN(LINKLINE); return LINK; }
<LINKLINE>{
  ("START"|"S"){WHITESPACE}*"="    { return START_NODE; }
  ("END"|"E"){WHITESPACE}*"="      { return END_NODE; }
  ("WORD"|"W"){WHITESPACE}*"="     { BEGIN(NWSTR); return WORD; }
  ("div"|"d"){WHITESPACE}*"="      { return DIV; }
  ("acoustic"|"a"){WHITESPACE}*"=" { return ACOUSTIC; }
  ("ngram"|"n"){WHITESPACE}*"="    { return NGRAM; }
  ("language"|"l"){WHITESPACE}*"=" { return LANGUAGE; }
  ("lmin"){WHITESPACE}*"="         { return LMIN; }
  ("lmout"){WHITESPACE}*"="        { return LMOUT; }
  ("posterior"|"p"){WHITESPACE}*"=" { return POSTERIOR; }
  ("error"|"e"){WHITESPACE}*"="    { return EXPECTED_ERROR; }
  ("correct"|"c"){WHITESPACE}*"="  { return CORRECT; }
  ("estimated_correct"|"c*"){WHITESPACE}*"="    { return ESTIMATED_CORRECT; }
  ("features"|"f"){WHITESPACE}*"=" { BEGIN(NWSTR); return FEATURES; }
  "x"{INTEGER}{WHITESPACE}*"="     { errno = 0;
                                    long n = strtol(yytext + 1,NULL,10); 
                                    if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                                      htk_error(yylloc, yyextra, "integer is out of range");
                                    }
                                    yylval->ival = n;
                         
                                    return XSCORE; 
                                  }
  ("output"|"O"){WHITESPACE}*"="    { BEGIN(NWSTR); return OUTPUT; }
  "#".*                            { /* yylloc->lines (1); ylloc->step (); return(COMMENT);*/ }
  {WHITESPACE}                     { /* yylloc->step(); */ }
  {NEWLINE}                        { END(); 
                                     // yylloc->lines (yyleng); yylloc->step (); 
                                     return ENDL;
                                   }
}


 
{NEWLINE}       { 
                  // yylloc->lines (yyleng); yylloc->step (); 
                  return ENDL;
                }

{WHITESPACE}    { /* yylloc->step(); *//*return(SPACE);*/ }

<NWSTR>{
  {NEWLINE}     { 
                  yylval->sval = htk_txtcpy("", false);
                  unput('\n'); 
                  END();
                  return STRING; 
                }

  {WHITESPACE}  { 
                  yylval->sval = htk_txtcpy("", false); 
                  END();
                  return STRING; 
                }

  {NOT_WSCHAR}* {
                  yylval->sval = htk_txtcpy(yytext, false); 
                  END();
                  return STRING; 
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
              return INT;
  		}
  
  {HEXADECIMAL}   {
  	          errno = 0;
  	          long n = strtol(yytext+2,NULL,16); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                htk_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
  	          return INT;
  		}
  
  {INTEGER}       { 
  	          errno = 0;
  	          long n = strtol(yytext,NULL,10); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                htk_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
              return INT;
  		}
  
  
  {FLOAT}         { 
  	          errno = 0;
  	          double n = strtod(yytext,NULL); 
              if (! (errno != ERANGE)) {
                htk_error(yylloc, yyextra, "float is out of range");
              }
  	          yylval->dval = n;
              return FLOAT;
  		}
}

<INITIAL,NODELINE,LINKLINE>{STRING} |
{CHARLIT}       { // Copy and remove "" and '' 
                  yylval->sval = htk_txtcpy(yytext, true); 
                  return STRING; 
                }


<*>{IDENT}{WHITESPACE}*"=" { 
                  yylval->sval = htk_txtcpy(yytext, false); 
                  BEGIN(NWSTR); 
                  return UNK_OPTION; 
                }

<*>.            { return ERR; }

%%

void HtkContext::init_scanner()
{
	yylex_init(&scanner);
	yyset_extra(this, scanner);
}

void HtkContext::destroy_scanner()
{
	yylex_destroy(scanner);
  scanner = 0;
}


