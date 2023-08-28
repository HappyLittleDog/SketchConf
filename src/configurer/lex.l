%option noyywrap
%option nounput
%option noinput
%{
#include <iostream>
#include <string>
#include "parser.tab.hpp"
#include "include/logger.hpp"

%}
%option yylineno
WhiteSpace              [ \t\n\r]*
LineComment             "//".*
NonZeroDigit            [1-9]
DecimalDigit            [0-9]
IDNonDigit              [A-Za-z_]
IDChar                  [0-9A-Za-z_]
%%
{WhiteSpace}                                            {;}
{LineComment}                                           {;}
"sketch"                                                {return SKETCH;}
"type"                                                  {return TYPE;}
"Pr"                                                    {return PROB;}
"err"                                                   {return ERROR;}
"CM"                                                    {yylval.type=0; return SKETCH_TYPE;}
"COUNT"                                                 {yylval.type=1; return SKETCH_TYPE;}
"TOWER"                                                 {yylval.type=2; return SKETCH_TYPE;}
"NITROCM"                                               {yylval.type=3; return SKETCH_TYPE;}
"NITROCOUNT"                                            {yylval.type=4; return SKETCH_TYPE;}
"RESULT"                                                {return RESULT;}
"ncols"                                                 {return NCOLS;}
"nrows"                                                 {return NROWS;}
"total memory"                                          {return MEM;}
"sample rate"                                           {return RATE;}
"\033[0m\033[1;35m["                                    {;}
"] \033[0m"                                             {;}
[:,<>=()]                                               {return yytext[0];}
"0"                                                     {yylval.val=0; return CONST;}
{NonZeroDigit}{DecimalDigit}*                           {yylval.val=std::strtod(yytext, nullptr); return CONST;}
{NonZeroDigit}{DecimalDigit}*\.{DecimalDigit}*          {yylval.val=std::strtod(yytext, nullptr); return CONST;}
0\.{DecimalDigit}*                                      {yylval.val=std::strtod(yytext, nullptr); return CONST;}
{IDNonDigit}{IDChar}*                                   {yylval.str_val=new std::string(yytext); return IDENT;}
.                                                       {LOG_ERROR("\"%s\"@[line #%d]: unrecognized token", yytext, yylineno);}
%%