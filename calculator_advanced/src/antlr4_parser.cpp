#include "../include/antlr4_parser.h"
#include "../grammar/MathLexer.h"
#include "../grammar/MathParser.h"
#include "../grammar/MathBaseVisitor.h"
#include <antlr4-runtime.h>
#include <iostream>
#include <stdexcept>

using namespace antlr4;

// Custom error listener implementation
class ANTLR4Parser::ErrorListener : public BaseErrorListener {
public:
    std::string lastError;
    
    void syntaxError(
        Recognizer *recognizer,
        Token *offendingSymbol,
        size_t line,
        size_t charPositionInLine,
        const std::string &msg,
        std::exception_ptr e) override
    {
        lastError = "Line " + std::to_string(line) + ":" + 
                    std::to_string(charPositionInLine) + " " + msg;
    }
};

// Expression visitor implementation
class ExpressionVisitor : public MathBaseVisitor {
private:
    ParsedEquation& equation_;
    
public:
    explicit ExpressionVisitor(ParsedEquation& eq) : equation_(eq) {}
    
    std::any visitNumber(MathParser::NumberContext *ctx) override {
        equation_.numbers.push_back(ctx->getText());
        return visitChildren(ctx);
    }
    
    std::any visitVariable(MathParser::VariableContext *ctx) override {
        equation_.variables.push_back(ctx->getText());
        return visitChildren(ctx);
    }
    
    std::any visitFunctionCall(MathParser::FunctionCallContext *ctx) override {
        equation_.functions.push_back(ctx->FUNCTION()->getText());
        return visitChildren(ctx);
    }
    
    std::any visitDerivativeExpr(MathParser::DerivativeExprContext *ctx) override {
        equation_.type = EquationType::DERIVATIVE;
        return visitChildren(ctx);
    }
    
    std::any visitIntegralExpr(MathParser::IntegralExprContext *ctx) override {
        equation_.type = EquationType::INTEGRAL;
        return visitChildren(ctx);
    }
    
    std::any visitSeriesExpr(MathParser::SeriesExprContext *ctx) override {
        equation_.type = EquationType::SERIES;
        return visitChildren(ctx);
    }
    
    std::any visitParametricExpr(MathParser::ParametricExprContext *ctx) override {
        equation_.type = EquationType::PARAMETRIC;
        return visitChildren(ctx);
    }
    
    std::any visitOdeExpr(MathParser::OdeExprContext *ctx) override {
        equation_.type = EquationType::ODE;
        return visitChildren(ctx);
    }
};

ANTLR4Parser::ANTLR4Parser() 
    : errorListener_(std::make_unique<ErrorListener>())
{
}

ANTLR4Parser::~ANTLR4Parser() = default;

ParsedEquation ANTLR4Parser::parse(const std::string& input) {
    ParsedEquation result;
    result.rawInput = input;
    result.isValid = false;
    
    try {
        // Create input stream
        ANTLRInputStream inputStream(input);
        
        // Create lexer
        MathLexer lexer(&inputStream);
        CommonTokenStream tokens(&lexer);
        
        // Create parser
        MathParser parser(&tokens);
        parser.removeErrorListeners();
        parser.addErrorListener(errorListener_.get());
        
        // Parse expression
        tree::ParseTree* tree = parser.start();
        
        // Check for syntax errors
        if (parser.getNumberOfSyntaxErrors() > 0) {
            result.errorMessage = errorListener_->lastError;
            return result;
        }
        
        // Visit tree to extract info
        ExpressionVisitor visitor(result);
        visitor.visit(tree);
        
        result.isValid = true;
        
    } catch (const std::exception& e) {
        result.errorMessage = std::string("Parse error: ") + e.what();
    }
    
    return result;
}

ParsedEquation ANTLR4Parser::parseUQFF(const std::string& input) {
    ParsedEquation result = parse(input);
    
    if (result.isValid) {
        result.isUQFF = true;
        
        // Check for UQFF-specific patterns
        if (input.find("F_U_Bi_i") != std::string::npos ||
            input.find("Um(") != std::string::npos ||
            input.find("g_MUGE") != std::string::npos) {
            result.type = EquationType::UQFF;
        }
    }
    
    return result;
}

bool ANTLR4Parser::validate(const std::string& input) {
    ParsedEquation result = parse(input);
    return result.isValid;
}

std::string ANTLR4Parser::getLastError() const {
    return errorListener_->lastError;
}
