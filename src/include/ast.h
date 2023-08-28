#pragma once
#include "defs.hpp"
#include "logger.hpp"
#include "configurer.h"
#include <string>
#include <deque>
#include <memory>

class BaseAst
{
public:
    virtual ~BaseAst() = default;
    virtual void Print(std::string indent = "") = 0;
};

class ASTDemands;
class ASTDemand;
class ASTConfigs;
class ASTConfig;
class ASTConstraints;
class ASTConstraint;

class ASTDemands : public BaseAst
{
public:
    std::deque<std::unique_ptr<ASTDemand> > items_;
    void Print(std::string indent = "") override;
    std::deque<Demand> Dump();
};

class ASTDemand : public BaseAst
{
public:
    std::string name_;
    SketchType type_;
    double rate_;
    std::unique_ptr<ASTConstraints> constraints_;
    void Print(std::string indent = "") override;
    Demand Dump();
};

class ASTConfigs : public BaseAst
{
public:
    std::deque<std::unique_ptr<ASTConfig> > items_;
    void Print(std::string indent = "") override;
    std::deque<Config> Dump();
};

class ASTConfig : public BaseAst
{
public:
    int type_;
    int nrows_;
    int ncols_;
    void Print(std::string indent = "") override;
    Config Dump();
};

class ASTConstraints : public BaseAst
{
public:
    std::deque<std::unique_ptr<ASTConstraint> > constraints_;
    void Print(std::string indent = "") override;
    std::deque<Constraint> Dump();
};

class ASTConstraint : public BaseAst
{
public:
    double err_;
    double prob_;
    void Print(std::string indent = "") override;
    Constraint Dump();
};