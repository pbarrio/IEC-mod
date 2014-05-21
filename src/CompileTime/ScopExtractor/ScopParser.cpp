/*
 * ScopParser.cpp: This file is part of the PolyOpt project
 *
 * PolyOpt: a Polyhedral Optimizer for the ROSE compiler
 *
 * Copyright (C) 2011 the Ohio State University
 *
 * This program can be redistributed and/or modified under the terms
 * of the license specified in the LICENSE.txt file at the root of the
 * project.
 *
 * Contact: P. Sadayappan <saday@cse.ohio-state.edu>
 *
 */
/**
 * @file: ScopParser.cpp
 * @author: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>, Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */


#include <rose.h>

#include <AstTraversal.h>
#include "ScopExtractor/ScopParser.hpp"
#include "ScopExtractor/SageNodeAnnotation.hpp"
#include "ScopExtractor/SageTools.hpp"


using namespace std;
using namespace SageInterface;


ScopParser::ScopParser()
{ }

ScopParser::~ScopParser()
{ }




/******************************************************************************/
/************************* High-level parser routines *************************/
/******************************************************************************/


/**
 * Check if a node respects the SCoP assumption.
 *
 *
 */
bool
ScopParser::isScopCompatible(SgNode* node,
			     std::set<SgNode*>& scopNodes,
			     std::set<SgNode*>& nonScopNodes)
{
  std::vector<SgNode*> readRefs;
  std::vector<SgNode*> writeRefs;

  if (isSgExprStatement(node))
    {
      // Polyhedral statement.
      SgExprStatement* stmt = isSgExprStatement(node);
      // Collect the references, and check one at most is a write reference.
      // Kind-a single-assignment.
      if (! collectReadWriteRefs(stmt, readRefs, writeRefs))
	return invalidNode(node, "Cannot collect refs for statement:");
      if (writeRefs.size() > 1)
	return invalidNode(node, "More than 1 write ref for statement:");

      // No side effect, see if all array accesses are affine.
      bool isScop = parseStatementExpression(stmt->get_expression());
      if (isScop)
	{
	  // Record the statement in the polyhedral statement list.
	  _polyStatements.insert(node);
	  // Add a specific annotation.
	  ScopStatementAnnotation* annot = new ScopStatementAnnotation();
	  annot->readRefs = readRefs;
	  annot->writeRefs = writeRefs;
	  node->setAttribute("ScopStmt", annot);
	  return validNode(node, "Scop OK:");
	}
      return invalidNode(node, "Statement is invalid: ");
    }

  if (isSgBasicBlock(node))
    {
      if (isBasicBlockNodeScopCompatible(isSgBasicBlock(node), scopNodes))
	return validNode(node, "Scop OK");
      return invalidNode(node, "Basic block is invalid: ");
    }

  if (isSgForStatement(node))
    {
      SgForStatement* forstmt = isSgForStatement(node);
      bool isScop = isForNodeScopCompatible(forstmt, scopNodes);

      if (isScop)
	{
	  //cout << forstmt->unparseToString() << endl;
	  // Record the for loop in the control stack.
	  updateControlStack(node);
	  return validNode(node, "Scop OK:");
	}
      return invalidNode(node, "FOR statement is invalid: ");
    }

  if (isSgIfStmt(node))
    {
      bool isScop = isIfNodeScopCompatible(isSgIfStmt(node), scopNodes);
      if (isScop)
	{
	  // Record the for loop in the control stack.
	  updateControlStack(node);
	  return validNode(node, "Scop OK: ");
	}
      return invalidNode(node, "IF statement is invalid: ");
    }


  if (isSgVariableDeclaration(node)) {
    // Be safe and lazy, return false.
    return invalidNode(node, "Variable declaration not supported in Scop");
  }

  // Unsupported node. Not a SCoP.
  return invalidNode(node, "Unsupported node type:");
}




/**
 * Check if a basic block is a scop.
 * - all statements in the BB are scops
 * - all affine expressions in the scop use read-only variables
 *
 */
bool
ScopParser::isBasicBlockNodeScopCompatible(SgBasicBlock* node,
					   std::set<SgNode*>& scopNodes) {

  // A BB is a scop if:
  // 1- all its children are scop.
  std::vector<SgNode*> children = node->get_traversalSuccessorContainer();
  std::vector<SgNode*>::const_iterator i;
  for (i = children.begin(); i != children.end(); ++i)
    {
      if (*i != NULL && scopNodes.find(*i) == scopNodes.end())
	break;
    }
  if (i != children.end())
    return false;

  // 2- all affine expressions are referencing read-only variables.
  std::vector<SgNode*> readRefsDummy;
  std::vector<SgNode*> writeRefsBB;
  if (! collectReadWriteRefs(isSgStatement(node),
			     readRefsDummy, writeRefsBB))
    return false;

  std::set<SgVariableSymbol*> writeRefsBBSymb =
    convertToSymbolSet(writeRefsBB);

  // For all children c, we have:
  // wrefs = writerefs(bb) - writerefs(c)
  // if wrefs inter affineref(c) not empty, then not a scop.
  // for (i = children.begin(); i != children.end(); ++i)
  //   {
  //     // Collect the r/w refs for the child.
  //     std::vector<SgNode*> writeRefsChild;
  //     if (! collectReadWriteRefs(isSgStatement(*i),
  // 				 readRefsDummy,
  // 				 writeRefsChild))
  // 	return false;

  //     std::set<SgVariableSymbol*> writeRefsChildSymb =
  // 	convertToSymbolSet(writeRefsChild);

  //     // Collect the refs in affine expr for the child.
  //     std::vector<SgNode*> affineRefsChild;
  //     collectAffineRefs(*i, affineRefsChild);

  //     std::set<SgVariableSymbol*> affineRefsChildSymb =
  // 	convertToSymbolSet(affineRefsChild);

  //     // Compute wrefs.
  //     std::set<SgVariableSymbol*> wrefs;
  //     set_difference(writeRefsBBSymb.begin(), writeRefsBBSymb.end(),
  // 		     writeRefsChildSymb.begin(), writeRefsChildSymb.end(),
  // 		     std::inserter(wrefs, wrefs.begin()));

  //     // Compute wrefs intersection with affineRef.
  //     std::set<SgVariableSymbol*> inter;
  //     set_intersection(wrefs.begin(), wrefs.end(),
  // 		       affineRefsChildSymb.begin(),
  // 		       affineRefsChildSymb.end(),
  // 		       std::inserter(inter, inter.begin()));

  //     if (inter.size())
  // 	return false;
  //   }

  // 3- No loop iterator is used outside its enclosing loop.
  std::vector<SgNode*> siblings;
  if (children.size())
    {
      i = children.begin();
      siblings.push_back(*i);
      ++i;
      for (; i != children.end(); ++i)
	{
	  if (! checkIteratorUsedOutsideLoop(siblings, *i))
	    return false;
	  else
	    siblings.push_back(*i);
	}
    }

  return true;
}


/**
 * Ensure a symbol written in the scop-to-be is not used in any affine
 * expression of the node to test, and conversely.
 * algo:
 *  if writtensymb(siblings U node) /\ affinesymb(siblings U node) != empty,
 *  then not a scop.
 *
 *
 */
bool ScopParser::checkAffineExpressionsReadOnlySymbols(std::vector<SgNode*>&
						       siblings,
						       SgNode* node)
{
  // Always true if nothing to process.
  if (siblings.size() == 0 || node == NULL)
    return true;

  // Collect all write references in the sibling list + node.
  std::vector<SgNode*>::const_iterator i;
  std::vector<SgNode*> readRefsDummy;
  std::vector<SgNode*> writeRefsNode;
  std::vector<SgNode*> writeRefsSiblings;
  std::vector<SgNode*> affineRefsNode;
  std::vector<SgNode*> affineRefsSiblings;
  for (i = siblings.begin(); i != siblings.end(); ++i)
    {
      if (! collectReadWriteRefs(isSgStatement(*i),
				 readRefsDummy, writeRefsNode))
	return false;
      writeRefsSiblings.insert(writeRefsSiblings.end(),
			       writeRefsNode.begin(), writeRefsNode.end());
      writeRefsNode.empty();
      readRefsDummy.empty();
      std::vector<SgNode*> affineRef;
      collectAffineRefs(*i, affineRef);
      affineRefsSiblings.insert (affineRefsSiblings.end(), affineRef.begin(),
				 affineRef.end());
    }
  std::set<SgVariableSymbol*> writeRefsSymb =
    convertToSymbolSet(writeRefsSiblings);
  std::set<SgVariableSymbol*> affineRefsSiblingsSymb =
    convertToSymbolSet(affineRefsSiblings);


  // Collect the r/w refs for the node.
  std::vector<SgNode*> writeRefsChild;
  if (! collectReadWriteRefs(isSgStatement(node),
			     readRefsDummy,
			     writeRefsChild))
    return false;

  std::set<SgVariableSymbol*> writeRefsChildSymb =
    convertToSymbolSet(writeRefsChild);
  // Collect the refs in affine expr for the node.
  std::vector<SgNode*> affineRefsChild;
  collectAffineRefs(node, affineRefsChild);
  std::set<SgVariableSymbol*> affineRefsChildSymb =
    convertToSymbolSet(affineRefsChild);


  // Prepare the union sets.
  std::set<SgVariableSymbol*> writtenSymbols = writeRefsSymb;
  writtenSymbols.insert (writeRefsChildSymb.begin(),
  			 writeRefsChildSymb.end());

  // Remove the loop iterators.
  std::vector<SgForStatement*> forloops =
    querySubTree<SgForStatement>(node, V_SgForStatement);
  for (i = siblings.begin(); i != siblings.end(); ++i)
    {
      std::vector<SgForStatement*> forloopss =
	querySubTree<SgForStatement>(*i, V_SgForStatement);
      forloops.insert(forloops.end(), forloopss.begin(), forloopss.end());
    }
  for (std::vector<SgForStatement*>::iterator i = forloops.begin();
       i != forloops.end(); ++i)
    {
      SgForStatement* fornode = isSgForStatement(*i);
      ROSE_ASSERT(fornode);
      ScopForAnnotation* annot =
	(ScopForAnnotation*)(fornode->getAttribute("ScopFor"));
      ROSE_ASSERT(annot);
      writtenSymbols.erase(annot->iterator);
    }
  std::set<SgVariableSymbol*> affinerefSymbols = affineRefsSiblingsSymb;
  affinerefSymbols.insert (affineRefsChildSymb.begin(),
			   affineRefsChildSymb.end());

  // Compute wrefs.
  // Compute wrefs /\ affineRef.
  std::set<SgVariableSymbol*> inter;
  set_intersection(writtenSymbols.begin(), writtenSymbols.end(),
		   affinerefSymbols.begin(),
		   affinerefSymbols.end(),
		   std::inserter(inter, inter.begin()));

  if (inter.size())
    return false;

  return true;
}



/**
 * Ensure a loop iterator value is not used outside the loop declaring
 * it.
 *
 *
 */
bool ScopParser::checkIteratorUsedOutsideLoop(std::vector<SgNode*>&
					      siblings,
					      SgNode* node)
{
  // Always true if nothing to process.
  if (siblings.size() == 0 || node == NULL)
    return true;

  // 1- Collect WS all loop iterator symbols in the sibling list.
  std::set<SgVariableSymbol*> loopIterators;
  std::vector<SgNode*>::const_iterator s;
  for (s = siblings.begin(); s != siblings.end(); ++s)
    {
      std::vector<SgNode*> loops =
	NodeQuery::querySubTree(*s, V_SgForStatement);
      std::vector<SgNode*>::const_iterator i;
      for (i = loops.begin(); i != loops.end(); ++i)
	{
	  SgForStatement* fornode = isSgForStatement(*i);
	  ROSE_ASSERT(fornode);
	  ScopForAnnotation* annot =
	    (ScopForAnnotation*)(fornode->getAttribute("ScopFor"));
	  ROSE_ASSERT(annot);
	  loopIterators.insert(annot->iterator);
	}
    }

  // 2- Collect WN all symbols read only in 'node'
  std::set<SgVariableSymbol*> readRefsSymb;
  collectReadOnlySymbols(isSgStatement(node), readRefsSymb);

  // 3- If WN /\ WS \ne \emptyset then a loop iterator value is used
  // outside its defining context.
  std::set<SgVariableSymbol*>::const_iterator i;
  for (i = readRefsSymb.begin(); i != readRefsSymb.end(); ++i)
    if (loopIterators.find(*i) != loopIterators.end())
      return false;

  return true;
}


/**
 * Check if a loop iterator value is used before the current scop.
 *
 */
bool ScopParser::checkIteratorUsedBeforeLoop(std::vector<SgNode*>& siblings,
					     SgNode* node)
{
  // Always true if nothing to process.
  if (siblings.size() == 0 || node == NULL)
    return true;

  // Do the converse operation to checkIteratorUsedOutsideLoop: look
  // for an iterator of 'node' being used in 'siblings'.

  std::vector<SgNode*>::const_iterator s;
  // 4- Collect WS all loop iterator symbols in the node.
  std::set<SgVariableSymbol*> nodeLoopIterators;
  std::vector<SgNode*> loops = NodeQuery::querySubTree(node, V_SgForStatement);
  std::vector<SgNode*>::const_iterator ii;
  for (ii = loops.begin(); ii != loops.end(); ++ii)
    {
      SgForStatement* fornode = isSgForStatement(*ii);
      ROSE_ASSERT(fornode);
      ScopForAnnotation* annot =
	(ScopForAnnotation*)(fornode->getAttribute("ScopFor"));
      ROSE_ASSERT(annot);
      nodeLoopIterators.insert(annot->iterator);
    }

  // 5- Collect WN all symbols read only in 'siblings'
  std::set<SgVariableSymbol*> siblingsReadRefsSymb;
  for (s = siblings.begin(); s != siblings.end(); ++s)
    collectReadOnlySymbols(isSgStatement(*s), siblingsReadRefsSymb);

  // 6- If WN /\ WS \ne \emptyset then a loop iterator value is used
  // outside its defining context.
  std::set<SgVariableSymbol*>::const_iterator i;
  for (i = siblingsReadRefsSymb.begin(); i != siblingsReadRefsSymb.end(); ++i)
    {
      if (nodeLoopIterators.find(*i) != nodeLoopIterators.end())
	return false;
    }

  return true;
}



/**
 * Check if an 'if' node is a scop.
 * - conditional is an affine expression
 * - no 'else' clause
 * - 'then' clause is a scop
 *
 */
bool
ScopParser::isIfNodeScopCompatible(SgIfStmt* ifnode,
				   std::set<SgNode*>& scopNodes)
{
  // Get the expression associated with the conditional.
  SgExprStatement* estmtcond = isSgExprStatement(ifnode->get_conditional());
  SgExpression* conditional = NULL;
  if (estmtcond != NULL)
    conditional = estmtcond->get_expression();

  // The 'then' clause must be a scop.
  if (scopNodes.find(ifnode->get_true_body()) == scopNodes.end())
    return false;
  else{
    //If there is an "else" clause, that must also be a scop
    if( ifnode->get_false_body() )
      if( scopNodes.find(ifnode->get_false_body()) == scopNodes.end() )
	return false;
  }

  // Collect read-only symbols dominated by the 'if' node.
  std::set<SgVariableSymbol*> symbols;
  collectReadOnlySymbols(ifnode, symbols);

  // Ensure the conditional is affine.
  // For the moment, no min/max, conjunctions of comparisons are ok.
  /// FIXME: Deal with min/max.
  bool res = parseExtendedAffineExpression(conditional);
    //isBoundExtendedAffineExpr(conditional, symbols);
  if (res)
    {
      ScopIfAnnotation* annot = new ScopIfAnnotation();
      annot->conditional = conditional;
      ifnode->setAttribute("ScopIf", annot);
    }

  return res;
}


/**
 * Check if a SgForNode respects the SCoP definition.
 * - loop has stride one (for the moment)
 * - loop body  is a scop
 * - loop iterator is not written in the body
 * - lb and ub of loop are affine expressions
 *
 */
bool
ScopParser::isForNodeScopCompatible(SgForStatement* fornode,
				    std::set<SgNode*>& scopNodes)
{
  // Retrieve the loop information.
  SgForInitStatement* lb = fornode->get_for_init_stmt();
  SgStatementPtrList initList = lb->get_init_stmt();
  SgStatement* ub = fornode->get_test();
  SgExpression* increment = fornode->get_increment();
  SgStatement* body = fornode->get_loop_body();

  // 1- Check the lb is of the form 'iterator = expression'.
  if (initList.size() != 1)
    return false;
  SgVariableSymbol* loopIterator = NULL;
  SgExpression* elb = NULL;
  
  SgStatement* init = *(initList.begin());
  // Allow C99 iterator declaration in for loop init.
  SgVariableDeclaration* vd = isSgVariableDeclaration(init);
  SgVariableSymbol* loopIteratorInDecl = NULL;
  if (vd)
    {
      return false;
    }
  else
    {
      parseAssignment(init, &loopIterator, &elb);
      if (loopIterator == NULL || elb == NULL)
	return false;
    }
    

  // 2- Loop has stride one.
  if (! parseSimpleIncrement(&loopIterator, increment))
    return false;

  // 3- Body is a scop.
  if (scopNodes.find(body) == scopNodes.end())
    return false;

  // 4- Bounds are affine expressions.

  // a) Get the ub expression.
  SgExprStatement* estmtub = isSgExprStatement(ub);
  SgExpression* eub = NULL;
  if (estmtub == NULL)
    return false;
  else{
    eub = estmtub->get_expression();
    if( ! (isSgLessThanOp(eub) && isSgVarRefExp(isSgLessThanOp(eub)->get_lhs_operand()) && (isSgVarRefExp(isSgLessThanOp(eub)->get_lhs_operand())->get_symbol() == loopIterator) ) )
      return false;
  }

  // b) Create the list of valid symbols for the expression.
  std::set<SgVariableSymbol*> symbols;
  // Side-effect analysis may fail. If so, we cannot determine which
  // symbols are read only, so by ROSE definition 'symbols' will be
  // empty.
  //collectReadWriteRefs(fornode, symbols);

  // Ensure the loop iterator is not written in the loop body.
  std::set<SgVariableSymbol*> bodySymbols;
  collectReadOnlySymbols(body, bodySymbols);
  if (bodySymbols.find(loopIterator) == bodySymbols.end())
    {
      // Visit the sub-tree, collecting all referenced variables.
      class ReferenceVisitor : public AstPrePostOrderTraversal {
      public:
	virtual void preOrderVisit(SgNode* node) {};
	virtual void postOrderVisit(SgNode* node) {
	  if (isSgVarRefExp(node) && ! isSgDotExp(node->get_parent()))
	    {
	      SgVarRefExp* vr = isSgVarRefExp(node);
	      SgVariableSymbol* symbol = vr->get_symbol();
	      symbols.insert(symbol);
	    }
	  else if (isSgDotExp(node))
	    {
	      SgVariableSymbol* symbol =
		SageTools::getFakeSymbol(isSgDotExp(node), fakeSymbolMap);
	      symbols.insert(symbol);
	    }
	}
	std::set<SgVariableSymbol*> symbols;
	std::map<std::string, SgVariableSymbol*> fakeSymbolMap;
      };

      ReferenceVisitor lookupReference;
      lookupReference.fakeSymbolMap = _fakeSymbolMap;
      lookupReference.traverse(body);
      // The loop iterator is written somewhere in the loop, or the
      // side-effect analysis cannot conclude it is read-only.
      if (lookupReference.symbols.find(loopIterator) !=
	  lookupReference.symbols.end())
	return false;
    }

  // c) Check the lb expression. 'max' are allowed but not 'min'.
  if (! parseExtendedAffineExpression(elb) )
    return false;

  // d) List of valid symbols now contains also the loop iterator for
  // the ub..
  //symbols.insert(loopIterator);

  // e) Check the ub expression. 'min' are allowed but not 'max'.
  if (! parseExtendedAffineExpression(eub))
    return false;

  ScopForAnnotation* annot = new ScopForAnnotation();
  annot->lb = elb;
  annot->ub = eub;
  annot->iterator = loopIterator;
  fornode->setAttribute("ScopFor", annot);

  return true;
}





/*****************************************************************************/
/*************************** Parsing functions *******************************/
/*****************************************************************************/





/**
 * Return true if the expression is an affine form of the symbols
 * referenced in the 'symbols' set.
 *
 */
bool
ScopParser::isAffineExpr(SgExpression* expr,
			 std::set<SgVariableSymbol*>& symbols)
{
  // 1- Ensure there are only standard binary operators in the expression.

  // Visit the sub-tree, collecting all referenced variables and
  // checking the type of nodes to match a candidate affine expression.
  class AffineVisitor : public AstPrePostOrderTraversal {
    public:
      virtual void preOrderVisit(SgNode* node) {};
    // Bottom-up parsing.
    virtual void postOrderVisit(SgNode* node) {
	if (! (isSgAddOp(node) ||
	       isSgIntegerDivideOp(node) ||
	       isSgMultiplyOp(node) ||
	       isSgSubtractOp(node) ||
	       isSgVarRefExp(node) ||
	       isSgDotExp(node) ||
	       isSgCastExp(node) ||
	       isSgMinusOp(node) ||
	       isSgUnaryAddOp(node) ||
	       SageTools::isIntegerTypeValue(isSgValueExp(node))))
	  hasMathBinopOnly = false;
	if (isSgDotExp(node))
	  {
	    SgVariableSymbol* symbol = parser->getFakeSymbol(isSgDotExp(node));
	    if (readSymbols.find(symbol) == readSymbols.end())
	      readSymbols.insert(symbol);
	  }
	else if (isSgVarRefExp(node) && ! isSgDotExp(node->get_parent()))
	  {
	    SgVarRefExp* vr = isSgVarRefExp(node);
	    SgVariableSymbol* symbol = vr->get_symbol();
	    if (readSymbols.find(symbol) == readSymbols.end())
	      readSymbols.insert(symbol);
	  }
      }
    std::set<SgVariableSymbol*> readSymbols;
    ScopParser* parser;
    bool hasMathBinopOnly;
    };

  AffineVisitor bv; bv.hasMathBinopOnly = true; bv.parser = this;
  bv.traverse(expr);

  if (bv.hasMathBinopOnly == false)
    return false;

  // 2- Ensure all referenced vars are in the allowed symbol list.
  std::set<SgVariableSymbol*>::const_iterator i;
  for (i = bv.readSymbols.begin(); i != bv.readSymbols.end(); ++i)
    if (symbols.find(*i) == symbols.end())
      return false;

  // 3- Check it is affine.
  return parseAffineExpression(expr, true);
}

/**
 * Return true if the expression is an affine form of the symbols
 * referenced in the 'symbols' set.
 *
 */
bool
ScopParser::isBoundExtendedAffineExpr(SgExpression* expr, std::set<SgVariableSymbol*>& symbols)
{
  // 1- Ensure there are only standard binary operators in the expression.

  // Visit the sub-tree, collecting all referenced variables and
  // checking the type of nodes to match a candidate affine expression.
  class AffineVisitor : public AstPrePostOrderTraversal {
    public:
      virtual void preOrderVisit(SgNode* node) {};
    // Bottom-up parsing.
    virtual void postOrderVisit(SgNode* node) {
      std::vector<SgFunctionSymbol*> call_stack;
      if (isSgArrowExp(node) ||
	  isSgArrowStarOp(node) ||
	  isSgAssignOp(node) ||
	  isSgCommaOpExp(node) || 
	  isSgCompoundAssignOp(node) ||
	  isSgConcatenationOp(node) ||
	  isSgDotExp(node) ||
	  isSgDotStarOp(node) ||
	  isSgIsNotOp(node) ||
	  isSgIsOp(node) ||
	  isSgMembershipOp(node) ||
	  isSgNonMembershipOp(node) ||
	  isSgPointerAssignOp(node) ||
	  isSgScopeOp(node) ||
	  isSgAddressOfOp(node) ||
	  isSgExpressionRoot(node) ||
	  isSgMinusMinusOp(node) ||
	  isSgPlusPlusOp(node) ||
	  isSgPointerDerefExp(node) ||
	  isSgThrowOp(node) ||
	  isSgUnaryAddOp(node) ||
	  ( isSgFunctionCallExp(node) && !SageTools::isMathFunc(isSgFunctionCallExp(node)->getAssociatedFunctionSymbol(),call_stack) )
	  )
	hasSideEffect = true;
      if (isSgDotExp(node)){
	SgVariableSymbol* symbol = parser->getFakeSymbol(isSgDotExp(node));
	if (readSymbols.find(symbol) == readSymbols.end())
	    readSymbols.insert(symbol);
      }
      else if (isSgVarRefExp(node) && ! isSgDotExp(node->get_parent())) {
	SgVarRefExp* vr = isSgVarRefExp(node);
	SgVariableSymbol* symbol = vr->get_symbol();
	if (readSymbols.find(symbol) == readSymbols.end())
	      readSymbols.insert(symbol);
      }
      else if( isSgPntrArrRefExp(node) ){
	SgExpression* lhs_exp = isSgPntrArrRefExp(node)->get_lhs_operand();
	if( !isSgVarRefExp(lhs_exp) )
	  indArray2D = true;
      }
    }
    std::set<SgVariableSymbol*> readSymbols;
    ScopParser* parser;
    bool hasSideEffect;
    bool indArray2D;
  };

  if( expr == NULL )
    return true;

  AffineVisitor bv; bv.hasSideEffect = false; bv.parser = this;
  bv.indArray2D = false;
  bv.traverse(expr);

  if (bv.hasSideEffect || bv.indArray2D)
    return false;

  // 2- Ensure all referenced vars are in the allowed symbol list.
  std::set<SgVariableSymbol*>::const_iterator i;
  for (i = bv.readSymbols.begin(); i != bv.readSymbols.end(); ++i)
    if (symbols.find(*i) == symbols.end())
      return false;

  // 3- Check it is affine.
  //return parseAffineExpression(expr, true);
  return true;
}


/**
 * Top-down parsing of an affine expression, returns true if it is
 * affine.
 * - var_ok = true -> expression can contain variable references.
 *
 */
bool
ScopParser::parseAffineExpression(SgExpression* expr, bool var_ok)
{
  SgBinaryOp* binop = isSgBinaryOp(expr);

  if (isSgAddOp(expr) || isSgSubtractOp(expr))
    {
      SgExpression* lhs = binop->get_lhs_operand();
      SgExpression* rhs = binop->get_rhs_operand();;
      bool lhs_val = parseAffineExpression(lhs, var_ok);
      bool rhs_val = parseAffineExpression(rhs, var_ok);
      return lhs_val && rhs_val;
    }
  if (isSgMultiplyOp(expr))
    {
      SgExpression* lhs = binop->get_lhs_operand();
      SgExpression* rhs = binop->get_rhs_operand();;
      if (NodeQuery::querySubTree(lhs, V_SgVarRefExp).size() != 0)
	{
	  bool lhs_val = parseAffineExpression(lhs, var_ok);
	  bool rhs_val = parseAffineExpression(rhs, false);
	  return lhs_val && rhs_val;
	}
      else if (NodeQuery::querySubTree(rhs, V_SgVarRefExp).size() != 0)
	{
	  bool lhs_val = parseAffineExpression(lhs, false);
	  bool rhs_val = parseAffineExpression(rhs, var_ok);
	  return lhs_val && rhs_val;
	}
      else
	{
	  bool lhs_val = parseAffineExpression(lhs, var_ok);
	  bool rhs_val = parseAffineExpression(rhs, var_ok);
	  return lhs_val && rhs_val;
	}
    }
  SgUnaryOp* unop = isSgUnaryOp(expr);
  if (isSgMinusOp(expr) || isSgUnaryOp(expr))
    return parseAffineExpression(unop->get_operand(), var_ok);

  /// FIXME: isIntegerDivideOp is not handled yet.
  if (SageTools::isIntegerTypeValue(isSgValueExp(expr)))
    return true;
  if (var_ok && getSymbolFromReference(expr))
    return true;

  return false;
}


/**
 * Top-down parsing of an extended affine expression, returns true if it is
 * affine.
 * - var_ok = true -> expression can contain variable references.
 *
 */

bool
ScopParser::parseExtendedAffineExpression(SgExpression* expr)
{
  SgBinaryOp* binop = isSgBinaryOp(expr);

  if( binop ){
    //Return false for side-effect expressions
    if( isSgArrowExp(binop) ||
	isSgArrowStarOp(binop) ||
	isSgAssignOp(binop) ||
	isSgCommaOpExp(binop) || 
	isSgCompoundAssignOp(binop) ||
	isSgConcatenationOp(binop) ||
	isSgDotExp(binop) ||
	isSgDotStarOp(binop) ||
	isSgIsNotOp(binop) ||
	isSgIsOp(binop) ||
	isSgMembershipOp(binop) ||
	isSgNonMembershipOp(binop) ||
	isSgPointerAssignOp(binop) ||
	isSgScopeOp(binop)
	)
      return false;
    // For Array access expressions, currently allowing only 1D arrays. Check if the access expression is itself "extended" affine.
    else if(isSgPntrArrRefExp(binop) ){
      SgExpression* lhs_exp = binop->get_lhs_operand();
      if( !isSgVarRefExp(lhs_exp) )
	return false;
      else
	return parseExtendedAffineExpression(binop->get_rhs_operand());
    }
    else{
      return parseExtendedAffineExpression(binop->get_lhs_operand()) && parseExtendedAffineExpression(binop->get_rhs_operand());
    }
  }
  SgUnaryOp* unaryop = isSgUnaryOp(expr);
  //Dissallow side-effect expression and pointer expressions
  if( unaryop ){
    if( isSgAddressOfOp(unaryop) ||
	isSgExpressionRoot(unaryop) ||
	isSgMinusMinusOp(unaryop) ||
	isSgPlusPlusOp(unaryop) ||
	isSgPointerDerefExp(unaryop) ||
	isSgThrowOp(unaryop) ||
	isSgUnaryAddOp(unaryop) 
	)
      return false;
    else
      return parseExtendedAffineExpression(unaryop->get_operand());
  }
  //Allow math functions
  SgFunctionCallExp* func_call_op = isSgFunctionCallExp(expr);
  if( func_call_op ){
    std::vector<SgFunctionSymbol*> call_stack;
    return SageTools::isMathFunc(func_call_op->getAssociatedFunctionSymbol(),call_stack);
  }
  //Of the remaining expression types allow only these expressions
  if ( isSgSizeOfOp(expr) ||
       isSgVarRefExp(expr) ||
       isSgValueExp(expr) )
    return true;
  else
    return false;  
}




/**
 * returns true if expr -> lhs '&&' rhs
 *
 *
 */
bool
ScopParser::parseConjunction(SgExpression* expr,
			     SgExpression** lhs,
			     SgExpression** rhs) {

  SgAndOp* op = isSgAndOp(expr);
  if (op)
    {
      *lhs = isSgExpression(op->get_lhs_operand());
      *rhs = isSgExpression(op->get_rhs_operand());
    }

  return true;
}


/**
 * returns true if stmt -> var '=' expr
 *
 *
 */
bool
ScopParser::parseAssignment(SgStatement* stmt, SgVariableSymbol** var,
			    SgExpression** expr)
{
  // Checking if node is of the form 'varname = <expression>'.
  SgExprStatement* assign = isSgExprStatement(stmt);
  if (! assign)
    return false;
  SgAssignOp* asop = isSgAssignOp(assign->get_expression());
  if (! asop)
    return false;
  *var = getSymbolFromReference(asop->get_lhs_operand());
  if (! *var)
    return false;
  *expr = isSgExpression(asop->get_rhs_operand());

  return true;
}


/**
 * returns true if:
 *       expr -> var++
 *       expr -> ++var
 *       expr -> var += 1
 *
 *
 */
bool ScopParser::parseSimpleIncrement(SgVariableSymbol** var,
				      SgExpression* expr)
{
  SgPlusPlusOp* ppop = isSgPlusPlusOp(expr);
  if (ppop)
    {
      SgVariableSymbol* symbol =
	getSymbolFromReference(ppop->get_operand_i());
      if (! symbol)
	return false;
      if (*var == NULL)
	*var = symbol;
      return symbol == *var;
    }
  SgPlusAssignOp* pop = isSgPlusAssignOp(expr);
  if (isSgPlusAssignOp(expr))
    {
      SgVariableSymbol* symbol =
	getSymbolFromReference(pop->get_lhs_operand());
      SgIntVal* val = isSgIntVal(pop->get_rhs_operand());
      if (! symbol)
	return false;
      if (*var == NULL)
	*var = symbol;
      return (symbol == *var && val && val->get_value() == 1);
    }

  return false;
}


/**
 * Generic expression parser. Returns true if the array access
 * functions are affine.
 *
 *
 *
 */
bool ScopParser::parseExpression(SgExpression* node) {
  if (isSgPntrArrRefExp(node))
    {
      SgPntrArrRefExp* arref;
      SgNode* n = node;
      // Grammar says we have 2 children, 1st is reference and 2nd is
      // offset. Hence offset must be an affine expression.

      bool isAffine = true;
      //SgPntrArrRefExp* last = arref;
      std::vector<bool> arrayAffineDim;
      while ((arref = isSgPntrArrRefExp(n)))
	{
	  bool isaff = parseAffineExpression(arref->get_rhs_operand(), true);
	  arrayAffineDim.push_back(isaff);
	  isAffine &= isaff;
	  n = arref->get_lhs_operand();
	}
      if (_nonaffineArrays && ! isAffine && n)
	{
	  SgVariableSymbol* symbol = getSymbolFromReference(n);
	  if (! symbol)
	    return invalidNode(node,
			       "Non affine access (no known approximattion):");

	  // Get the previous entry for the array, if any, in the list
	  // of non-affine arrays.
	  std::set<std::pair<SgVariableSymbol*, std::vector<bool> > >::iterator
	    i;
	  for (i = _arraysAsScalars.begin(); i != _arraysAsScalars.end(); ++i)
	    if (i->first == symbol)
	      break;
	  if (i != _arraysAsScalars.end())
	    {
	      for (int j = 0; j < arrayAffineDim.size(); ++j)
		arrayAffineDim[j] = arrayAffineDim[j] && i->second[j];
	      _arraysAsScalars.erase(i);
	    }
	  std::pair<SgVariableSymbol*, std::vector<bool> >
	    elt(symbol, arrayAffineDim);
	  _arraysAsScalars.insert(elt);
	  return validNode(node, "Non affine reference approximated: ");
	}
      if (isAffine)
	return validNode(node, "Affine array access: ");
      else
	return invalidNode(node,
			   "Non affine access (consider approximate SC?):");
    }
  if (isSgFunctionCallExp(node))
    {
      // Function has no side-effect, unless we wouldn't be there.
      // Process its arguments.
      SgFunctionCallExp* fc = isSgFunctionCallExp(node);
      SgExprListExp* args = fc->get_args();
      bool val = true;
      if (args)
	{
	  SgExpressionPtrList exprs = args->get_expressions();
	  SgExpressionPtrList::const_iterator i;
	  for (i = exprs.begin(); i != exprs.end(); ++i)
	    val = val && parseExpression(*i);

	}
      return val;
    }
  if (isSgBinaryOp(node))
    {
      SgBinaryOp* binop = isSgBinaryOp(node);
      SgExpression* lhs = binop->get_lhs_operand();
      SgExpression* rhs = binop->get_rhs_operand();;
      bool lhs_val = parseExpression(lhs);
      bool rhs_val = parseExpression(rhs);
      return lhs_val && rhs_val;
    }
  if (isSgUnaryOp(node))
    return parseExpression(isSgUnaryOp(node)->get_operand());

  // The rest is always ok, we have already checked the expression is
  // in single assignment form.
  return true;
}


/**
 * Generic expression parser. Returns true if the array access
 * functions are affine.
 *
 *
 *
 */
bool ScopParser::parseStatementExpression(SgExpression* node) {
  if (isSgPntrArrRefExp(node))
    {
      SgPntrArrRefExp* array_access_exp = isSgPntrArrRefExp(node);
      SgExpression* array_var = array_access_exp->get_lhs_operand();
      SgExpression* access_exp = array_access_exp->get_rhs_operand();
      // Grammar says we have 2 children, 1st is reference and 2nd is
      // offset. For inner dimensions the expression must be affine.
      // For outer-most dimension it can be non-affine
      
      bool isAffine = true;
      while( isSgPntrArrRefExp(array_var) ){
	isAffine &= parseAffineExpression(access_exp,true);
	array_access_exp = isSgPntrArrRefExp(array_var);
	array_var = array_access_exp->get_lhs_operand();
	access_exp = array_access_exp->get_rhs_operand();
      }
      if( isAffine ){
	bool isEAffine = parseExtendedAffineExpression(access_exp);
	if( isEAffine )
	  return validNode(node,"Outer_most dimension is extended affine");
	else
	  return invalidNode(node,"Outer most dimension is not extended affine");
      }
      else{
	return invalidNode(node,"Inner Array dimension is not affine");
      }
    }
  if (isSgFunctionCallExp(node))
    {
      // Function has no side-effect, unless we wouldn't be there.
      // Process its arguments.
      SgFunctionCallExp* fc = isSgFunctionCallExp(node);
      SgExprListExp* args = fc->get_args();
      bool val = true;
      if (args)
	{
	  SgExpressionPtrList exprs = args->get_expressions();
	  SgExpressionPtrList::const_iterator i;
	  for (i = exprs.begin(); i != exprs.end(); ++i)
	    val = val && parseStatementExpression(*i);

	}
      return val;
    }
  if (isSgBinaryOp(node))
    {
      SgBinaryOp* binop = isSgBinaryOp(node);
      SgExpression* lhs = binop->get_lhs_operand();
      SgExpression* rhs = binop->get_rhs_operand();;
      bool lhs_val = parseStatementExpression(lhs);
      bool rhs_val = parseStatementExpression(rhs);
      return lhs_val && rhs_val;
    }
  if (isSgUnaryOp(node))
    return parseStatementExpression(isSgUnaryOp(node)->get_operand());

  // The rest is always ok, we have already checked the expression is
  // in single assignment form.
  return true;
}



/**
 * Collect all SgVarRefExp nodes in the subtree, such that a reference
 * is being used in an affine expression (array access, loop bound or
 * conditional).
 *
 *
 */
bool
ScopParser::collectAffineRefs(SgNode* node,
			      std::vector<SgNode*>& refs) {
  std::vector<SgNode*> arrays =
    NodeQuery::querySubTree(node, V_SgPntrArrRefExp);

  std::vector<SgNode*> conditionals =
    NodeQuery::querySubTree(node, V_SgIfStmt);

  std::vector<SgNode*> loops =
    NodeQuery::querySubTree(node, V_SgForStatement);

  std::vector<SgNode*>::const_iterator i;
  std::vector<SgNode*> readRefs;
  std::vector<SgNode*> writeRefs;

  // Internal helper:
  // Get all references used in a given expression
  class GetVarRefsUsedInExpressionVisitor: public AstSimpleProcessing
  {
  public:

    virtual void visit(SgNode* n) {
      if (isSgVarRefExp(n))
	refs.insert(n);
    }
    std::set<SgNode*> refs;
  };


  // Collect read references for all array references.
  for (i = arrays.begin(); i != arrays.end(); ++i)
    {
      ROSE_ASSERT(*i);
      SgPntrArrRefExp* arref = isSgPntrArrRefExp(*i);
      ROSE_ASSERT(arref);
      GetVarRefsUsedInExpressionVisitor vis;
      vis.traverse(arref->get_rhs_operand(), preorder);

      std::set<SgNode*>::const_iterator j;
      for (j = vis.refs.begin(); j != vis.refs.end(); ++j)
	refs.push_back(*j);
    }

  // Collect read references for all conditionals.
  for (i = conditionals.begin(); i != conditionals.end(); ++i)
    {
      readRefs.empty();
      writeRefs.empty();
      ROSE_ASSERT(*i);
      SgIfStmt* ifnode = isSgIfStmt(*i);
      if (! collectReadWriteRefs(ifnode->get_conditional(),
				 readRefs, writeRefs))
	return false;
      ROSE_ASSERT(writeRefs.size() == 0);
      refs.insert(refs.end(), readRefs.begin(), readRefs.end());
    }

  // Collect read references for all loop bounds.
  for (i = loops.begin(); i != loops.end(); ++i)
    {
      readRefs.empty();
      writeRefs.empty();
      ROSE_ASSERT(*i != NULL);
      SgForStatement* fornode = isSgForStatement(*i);
      ROSE_ASSERT(fornode != NULL);
      if (! collectReadWriteRefs(fornode->get_for_init_stmt(),
				 readRefs, writeRefs))
	return false;
      refs.insert(refs.end(), readRefs.begin(), readRefs.end());

      readRefs.empty();
      writeRefs.empty();
      if (! collectReadWriteRefs(fornode->get_test(),
				 readRefs, writeRefs))
	return false;
      refs.insert(refs.end(), readRefs.begin(), readRefs.end());
    }

  return true;
}



/*****************************************************************************/
/************************** Interface to SageTools ***************************/
/*****************************************************************************/


/**
 * Interface to SageTools::invalidNode, using SageParser::_verboseLevel.
 *
 */
bool
ScopParser::invalidNode(SgNode* node, const char* message)
{
  //return SageTools::invalidNode(node, message);
  return false;
}

/**
 * Interface to SageTools::validNode, using SageParser::_verboseLevel.
 *
 */
bool
ScopParser::validNode(SgNode* node, const char* message)
{
  //return SageTools::validNode(node, message);
  return true;
}


/**
 * Interface to SageTools::convertToSymbolSet, using
 * SageTools::_fakeSymbolMap.
 *
 *
 */
std::set<SgVariableSymbol*>
ScopParser::convertToSymbolSet(std::vector<SgNode*>& refs)
{
  return SageTools::convertToSymbolSet(refs, _fakeSymbolMap);
}


/**
 * Interface to SageTools::collectReadOnlySymbols, using
 * SageTools::_fakeSymbolMap.
 *
 */
bool
ScopParser::collectReadOnlySymbols(SgStatement* node,
				   std::set<SgVariableSymbol*>& symbols)
{
  return SageTools::collectReadOnlySymbols(node, symbols, _fakeSymbolMap,
					   _allowMathFunc);
}

/**
 * Interface to SageTools::collectReadWriteRefs, using
 * SageTools::_fakeSymbolMap.
 *
 */
bool
ScopParser::collectReadWriteRefs(SgStatement* node,
				std::vector<SgNode*>& readRefs,
				std::vector<SgNode*>& writeRefs)
{
  return SageTools::collectReadWriteRefs(node, readRefs, writeRefs,
					 _allowMathFunc);
}

/**
 * Interface to SageTools::getFakeSymbol, using
 * SageParser::_fakeSymbolMap.
 *
 */
SgVariableSymbol*
ScopParser::getFakeSymbol(SgDotExp* node)
{
  return SageTools::getFakeSymbol(node, _fakeSymbolMap);
}

/**
 * Interface to SageTools::getSymbolFromReference, using
 * SageParser::_fakeSymbolMap.
 *
 */
SgVariableSymbol*
ScopParser::getSymbolFromReference(SgNode* node)
{
  return SageTools::getSymbolFromReference(node, _fakeSymbolMap);
}




/*****************************************************************************/
/*************************** Array delinearization ***************************/
/*****************************************************************************/

/**
 * Inspect all references to an array that has been marked as having a
 * non-affine access function. Check if the references represent a
 * linearization of a multidimensional access function.
 *
 */
void
ScopParser::checkLinearizedArrayAccesses(SgNode* root)
{
  // Iterate on all arrays with non-affine array access.
  std::set<std::pair<SgVariableSymbol*, std::vector<bool> > >::const_iterator i;
  for (i = _arraysAsScalars.begin(); i != _arraysAsScalars.end(); ++i)
    {
      // 1- Collect all references to this array.

      // 2- Ensure all references are to a 1-d array.

      // 3- Ensure all references are a polynomial form of the same
      // parameters, for all references.

      // 3.1- Collect all read-only symbols in the scop.
      // 3.2- Collect all symbols used in the array reference.
      // 3.3- The reference use only read-only symbols and
      // surrounding loop iterator. Required since we did no check on
      // the reference as it is non-affine.
      // 3.4- Build the polynomial as the form of the parameters.


    }

}




/*****************************************************************************/
/*************************** Control stack handling **************************/
/*****************************************************************************/

/**
 * Reset the parser state.
 *
 */
void
ScopParser::emptyParser() {
  _controlStack.empty();
  _polyStatements.empty();
  _controlMap.empty();
  _arraysAsScalars.empty();
  _fakeSymbolMap.empty();
  _linearizedArraysAsArrays.empty();
}


/**
 * Update the control stack of each polyhedral statements with the
 * new control.
 *
 */
bool
ScopParser::updateControlStack(SgNode* control) {
  std::set<SgNode*>::const_iterator i;
  for (i = _polyStatements.begin(); i != _polyStatements.end(); ++i)
    {
      if (SageTools::isDominatedBy(*i, control))
	{
	  std::map<SgNode*, std::vector<SgNode*> >::iterator elt =
	    _controlMap.find(*i);
	  if (elt != _controlMap.end())
	    elt->second.push_back(control);
	  else
	    {
	      std::vector<SgNode*> controlList;
	      controlList.push_back(control);
	      std::pair<SgNode*, std::vector<SgNode*> > newelt(*i, controlList);
	      _controlMap.insert(newelt);
	    }
	}
    }

  return true;
}


/**
 * Finalize control map, by creating map entries for polyhedral
 * statements without any control.
 *
 */
void
ScopParser::finalizeControlMap() {
  std::set<SgNode*>::const_iterator i;
  for (i = _polyStatements.begin(); i != _polyStatements.end(); ++i)
    {
      std::map<SgNode*, std::vector<SgNode*> >::iterator j;
      j = _controlMap.find(*i);
      if (j == _controlMap.end())
	{
	  std::vector<SgNode*> control;
	  std::pair<SgNode*, std::vector<SgNode*> > newelt(*i, control);
	  _controlMap.insert(newelt);
	}
    }
}


/*****************************************************************************/
/****************************** getter / setter ******************************/
/*****************************************************************************/


/**
 * Accessor.
 *
 */
std::map<SgNode*, std::vector<SgNode*> >
ScopParser::getControlMap() {
  return _controlMap;
}

/**
 * Accessor.
 *
 */
std::set<std::pair<SgVariableSymbol*, std::vector<bool> > >
ScopParser::getArraysAsScalars() {
  return _arraysAsScalars;
}


/**
 * Accessor.
 *
 */
std::map<std::string, SgVariableSymbol*>
ScopParser::getFakeSymbolMap()
{
  return _fakeSymbolMap;
}

/**
 * Accessor.
 *
 */
std::set<SgVariableSymbol*>
ScopParser::getLinearizedArraysAsArrays()
{
  return _linearizedArraysAsArrays;
}


/**
 * Setter.
 *
 */
void
ScopParser::setVerboseLevel(int level)
{
  _verboseLevel = level;
}

/**
 * Setter.
 *
 */
void
ScopParser::setNonaffineArrays(bool val)
{
  _nonaffineArrays = val;
}


/**
 * Setter.
 *
 */
void
ScopParser::setAllowMathFunc(bool val)
{
  _allowMathFunc = val;
}
