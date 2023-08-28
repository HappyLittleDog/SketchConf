#include "include/ast.h"
#include <string>
#include <iostream>
using namespace std;

void ASTDemands::Print(string indent)
{
    for (int i=0;i<items_.size();i++)
        items_[i]->Print(indent);
}

deque<Demand> ASTDemands::Dump()
{
    deque<Demand> nt;
    for (int i=0;i<items_.size();i++)
        nt.push_back(items_[i]->Dump());
    return nt;
}

void ASTDemand::Print(string indent)
{
    cout<<indent<<"Demand {\n";
    auto curindent=indent+"|\t";
    cout<<curindent<<"Name: "<<name_<<endl;
    cout<<curindent<<"Type: "<< static_cast<int>(type_)<<endl;
    cout<<curindent<<"Sample rate: "<<rate_<<endl;
    constraints_->Print(curindent);
    cout<<indent<<"}"<<endl;
}

Demand ASTDemand::Dump()
{
    return {name_, type_, constraints_->Dump(), rate_};
}

void ASTConfigs::Print(string indent)
{
    for (int i=0;i<items_.size();i++)
        items_[i]->Print(indent);
}

deque<Config> ASTConfigs::Dump()
{
    deque<Config> nt;
    for (int i=0;i<items_.size();i++)
        nt.push_back(items_[i]->Dump());
    return nt;
}

void ASTConfig::Print(string indent)
{
    cout<<indent<<"Config {\n";
    auto curindent=indent+"|\t";
    cout<<curindent<<"Type: "<< type_ <<endl;
    cout<<curindent<<"#rows: "<< nrows_ <<endl;
    cout<<curindent<<"#cols: "<< ncols_ <<endl;
    cout<<indent<<"}"<<endl;
}

Config ASTConfig::Dump()
{
    return {type_, nrows_, ncols_};
}

void ASTConstraints::Print(string indent)
{
    for (int i=0;i<constraints_.size();i++)
    {
        constraints_[i]->Print(indent);
    }
}

deque<Constraint> ASTConstraints::Dump()
{
    deque<Constraint> nt;
    for (int i=0;i<constraints_.size();i++)
        nt.push_back(constraints_[i]->Dump());
    return nt;
}

void ASTConstraint::Print(string indent)
{
    cout<<indent<<"Constraint: P(ERROR > "<<err_<<") < "<<prob_<<endl;
}

Constraint ASTConstraint::Dump()
{
    return {err_, prob_};
}