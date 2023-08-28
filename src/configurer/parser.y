%code requires {
  #include <memory>
  #include <string>
  #include "include/ast.h"
  #include "include/logger.hpp"
}
%{
#include <iostream>
#include <memory>
#include <string>
#include "include/ast.h"
#include "include/logger.hpp"

int yylex();
void yyerror(std::unique_ptr<BaseAst> &ast, const char *s);
%}
%parse-param { std::unique_ptr<BaseAst> &ast }
%union {
    int type;
    double val;
    std::string* str_val;
    BaseAst* ast_val;
}
%token SKETCH TYPE PROB ERROR NCOLS NROWS MEM RESULT RATE
%token <type> SKETCH_TYPE
%token <str_val> IDENT
%token <val> CONST
%type <ast_val> Demands Demand Configs Config Constraints Constraint
%%
CompUnits
    :   Demands
        {
            ast=std::unique_ptr<BaseAst>($1);
        }
    |   Configs
        {
            ast=std::unique_ptr<BaseAst>($1);
        }
    ;

Demands
    :   Demand
        {
            auto cur=new ASTDemands();
            cur->items_.push_back(std::unique_ptr<ASTDemand>(reinterpret_cast<ASTDemand*>($1)));
            $$=cur;
        }
    |   Demand Demands
        {
            auto cur=reinterpret_cast<ASTDemands*>($2);
            cur->items_.push_front(std::unique_ptr<ASTDemand>(reinterpret_cast<ASTDemand*>($1)));
            $$=cur;
        }
    ;

Demand
    :   SKETCH IDENT ':' TYPE '=' SKETCH_TYPE  Constraints
        {
            auto cur=new ASTDemand();
            cur->name_ = *($2);
            cur->type_ = static_cast<SketchType>($6);
            cur->rate_ = -1;
            cur->constraints_ = std::unique_ptr<ASTConstraints>(reinterpret_cast<ASTConstraints*>($7));
            $$=cur;
        }
    |   SKETCH IDENT ':' TYPE '=' SKETCH_TYPE ',' RATE '=' CONST  Constraints
        {
            auto cur=new ASTDemand();
            cur->name_ = *($2);
            cur->type_ = static_cast<SketchType>($6);
            cur->rate_ = $10;
            cur->constraints_ = std::unique_ptr<ASTConstraints>(reinterpret_cast<ASTConstraints*>($11));
            $$=cur;
        }
    ;

Configs
    :   Config
        {
            auto cur=new ASTConfigs();
            cur->items_.push_back(std::unique_ptr<ASTConfig>(reinterpret_cast<ASTConfig*>($1)));
            $$=cur;
        }
    |   Config Configs
        {
            auto cur=reinterpret_cast<ASTConfigs*>($2);
            cur->items_.push_front(std::unique_ptr<ASTConfig>(reinterpret_cast<ASTConfig*>($1)));
            $$=cur;
        }
    ;

Config
    :   RESULT TYPE '=' CONST ',' NCOLS '=' CONST ',' NROWS '=' CONST ',' MEM '=' CONST
        {
            auto cur=new ASTConfig();
            cur->type_ = $4;
            cur->nrows_ = $12;
            cur->ncols_ = $8;
            $$=cur;
        }
    |   RESULT TYPE '=' CONST ',' NROWS '=' CONST ',' NCOLS '=' CONST ',' MEM '=' CONST
        {
            auto cur=new ASTConfig();
            cur->type_ = $4;
            cur->nrows_ = $8;
            cur->ncols_ = $12;
            $$=cur;
        }
    |   RESULT TYPE '=' CONST ',' NCOLS '=' CONST ',' NROWS '=' CONST
        {
            auto cur=new ASTConfig();
            cur->type_ = $4;
            cur->nrows_ = $12;
            cur->ncols_ = $8;
            $$=cur;
        }
    |   RESULT TYPE '=' CONST ',' NROWS '=' CONST ',' NCOLS '=' CONST
        {
            auto cur=new ASTConfig();
            cur->type_ = $4;
            cur->nrows_ = $8;
            cur->ncols_ = $12;
            $$=cur;
        }
    ;

Constraints
    :
        {
            auto cur=new ASTConstraints();
            $$=cur;
        }
    |   ',' Constraint Constraints
        {
            auto cur=reinterpret_cast<ASTConstraints*>($3);
            cur->constraints_.push_front(std::unique_ptr<ASTConstraint>(reinterpret_cast<ASTConstraint*>($2)));
            $$=cur;
        }
    ;

Constraint
    :   PROB '(' ERROR '>' CONST ')' '<' CONST
        {
            auto cur=new ASTConstraint();
            cur->err_=$5;
            cur->prob_=$8;
            $$=cur;
        }
    ;
%%
void yyerror(std::unique_ptr<BaseAst> &ast, const char *s)
{
    LOG_ERROR("%s", s);
}
