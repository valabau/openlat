/* based on http://www.phpcompiler.org/articles/reentrantparser.html */
%option reentrant
%option prefix="openlat_thot_"
%option bison-bridge
%option bison-locations
%option noyywrap
%option yylineno
%option batch 
%option debug 
%option case-insensitive

%{

#include "thot.h"

using namespace openlat;
using namespace openlat::thot;


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
                thot_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
              return parser::token::INT;
  		}
  
  {HEXADECIMAL}   {
  	          errno = 0;
  	          long n = strtol(yytext+2,NULL,16); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                thot_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
  	          return parser::token::INT;
  		}
  
  {INTEGER}       { 
  	          errno = 0;
  	          long n = strtol(yytext,NULL,10); 
              if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE)) {
                thot_error(yylloc, yyextra, "integer is out of range");
              }
  	          yylval->ival = n;
              return parser::token::INT;
  		}
  
  
  {FLOAT}         { 
  	          errno = 0;
  	          double n = strtod(yytext,NULL); 
              if (! (errno != ERANGE)) {
                thot_error(yylloc, yyextra, "float is out of range");
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


<*>.            { thot_error(yylloc, yyextra, "invalid character"); exit(-1); }

%%

void openlat::thot::ThotContext::init_scanner()
{
	yylex_init(&scanner);
	yyset_extra(this, scanner);
}

void openlat::thot::ThotContext::destroy_scanner()
{
	yylex_destroy(scanner);
  	scanner = 0;
}

void openlat::thot::thot_error(YYLTYPE* locp, openlat::thot::ThotContext*, const char* err) 
{
	openlat::parser_error<YYLTYPE>(*locp, string(err));
}

void openlat::thot::parser::error(const YYLTYPE &locp, const string &err) 
{
	openlat::parser_error<YYLTYPE>(locp, err);
}