/**
 * @file Z3Tools.h
 * @author Vincent Penelle (vincent.penelle@u-bordeaux.fr)
 * @brief This file contains functions to easily and transparently manipulate the Z3 SMT solver as a SAT solver without delving too much in the documentation.
 *        Note that most functions ask for a Z3_context. This is needed by Z3, but in this case, simply create a context at the beginning and pass it as an argument
 *        every time needed.
 * @version 1
 * @date 2019-08-01
 * 
 * @copyright Creative Commons
 * 
 */

#ifndef COCA_Z3TOOLS_H_
#define COCA_Z3TOOLS_H_

#include <z3.h>
#include <stdbool.h>

/**
 * @brief Creates a basic Z3 context with basic config (sufficient for this project). Must be freed at end of program with Z3_del_context.
 * 
 * @return Z3_context The created context
 */
Z3_context makeContext(void);

/**
 * @brief Creates a formula containing a single variable whose name is given in parameter. Example mk_bool_var(ctx,"toto") will create the formula «toto». Each call with
 *        same name will produce the same formula (so it can be used to have the same variable in different formulae.)
 * 
 * @param ctx The context of the solver.
 * @param name The name the variable.
 * @return Z3_ast The formula consisting in the variable.
 */
Z3_ast mk_bool_var(Z3_context ctx, const char * name);

/**
 * @brief Tells if a formula is satisfiable, unsatisfiable, or cannot be decided.
 * 
 * @param ctx The context of the solver.
 * @param formula The formula to check.
 * @return Z3_lbool Z3_L_FALSE if @p formula is unsatisfiable, Z3_L_TRUE if @p formula is satisfiable and Z3_L_UNDEF if the solver cannot decide if @p formula is
 *         satisfiable or not.
 */
Z3_lbool isFormulaSat(Z3_context ctx, Z3_ast formula);

/**
 * @brief Returns an assignment of variables satisfying the formula if it is satisfiable. Exits the program if the formula is unsatisfiable.
 * 
 * @param ctx The context of the solver.
 * @param formula The formula to get a model from.
 * @return Z3_model An assignment of variable satisfying @p formula.
 */
Z3_model getModelFromSatFormula(Z3_context ctx, Z3_ast formula);

/**
 * @brief Checks if a formula is satisfiable, unsatisfiable, or cannot be decided. If it is decidable, puts a model in the formula in model.
 * 
 * @param ctx The context of the solver.
 * @param formula The formula to check.
 * @param model A pointer towards a model. Will contain a model of @p formula if it is satisfiable (otherwise, will not be modified).
 * @return Z3_lbool Z3_L_FALSE if @p formula is unsatisfiable, Z3_L_TRUE if @p formula is satisfiable and Z3_L_UNDEF if the solver cannot decide if @p formula is satisfiable or not.
 */
Z3_lbool solveFormula(Z3_context ctx, Z3_ast formula, Z3_model* model);

/**
 * @brief Returns the truth value of the formula @p variable in the variable assignment @p model. Very usefull if @p variable is a formula containing a single variable.
 * 
 * @param ctx The context of the solver.
 * @param model A variable assignment.
 * @param variable A formula of which we want the truth value.
 * @return true if @p variable is true in @p model
 * @return false otherwise.
 * @pre @p model must be a valid model.
 */
bool valueOfVarInModel(Z3_context ctx, Z3_model model, Z3_ast variable);

#endif