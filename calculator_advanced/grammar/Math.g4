grammar Math;

// Parser rules
start: expression EOF ;

expression: expr (';' expr)* ;

expr: NUMBER                                 # Number
    | VARIABLE                               # Variable
    | expr op=('*'|'/') expr                 # Mul
    | expr op=('+'|'-') expr                 # Add
    | expr '^' expr                          # Pow
    | '(' expr ')'                           # Paren
    | FUNCTION '(' expr (',' expr)* ')'      # FunctionCall
    | derivative                             # DerivativeExpr
    | integral                               # IntegralExpr
    | summation                              # SummationExpr
    | product                                # ProductExpr
    | series                                 # SeriesExpr
    | parametric                             # ParametricExpr
    | ode                                    # OdeExpr
    | expr '=' expr                          # Equation
    ;

derivative: ('d/d' VARIABLE | '∂/∂' VARIABLE)+ expr ;

integral: '∫' limits? expr 'd' VARIABLE ;

summation: '∑' limits? expr ;

product: '∏' limits? expr ;

series: 'series' '(' expr ',' VARIABLE ',' NUMBER ')' ;

parametric: VARIABLE '(' VARIABLE ')' '=' expr ;

ode: 'd' VARIABLE '/d' VARIABLE '=' expr ;

limits: '(' expr ',' expr ')' ;

// Lexer rules
NUMBER: [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;

VARIABLE: [a-zA-Z] [a-zA-Z0-9_]* ;

FUNCTION: 'sin' | 'cos' | 'tan' | 'exp' | 'log' | 'sqrt' | 'abs' ;

WS: [ \t\r\n]+ -> skip ;

COMMENT: '//' ~[\r\n]* -> skip ;
