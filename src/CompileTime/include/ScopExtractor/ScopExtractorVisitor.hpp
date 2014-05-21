/*
 * ScopExtractorVistor.hpp: This file is part of the PolyOpt project
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
 * @file: ScopExtractorVistor.hpp
 * @author: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>, Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */

#ifndef PLUTO_ROSE_SCOP_EXTRACTOR_VISITOR_HPP
# define PLUTO_ROSE_SCOP_EXTRACTOR_VISITOR_HPP

# include <AstTraversal.h>
# include "ScopExtractor/ScopParser.hpp"
// # include <sideEffect.h>


class ScopExtractorVisitor : public AstPrePostOrderTraversal {
public:
  // Special constructor to store the root node of the sub-tree.
  ScopExtractorVisitor(SgNode* root);

  virtual void preOrderVisit(SgNode* node) {}

  // Postfix algorithm for SCoP detection.
  virtual void postOrderVisit(SgNode* node);

  void createScopFromNodes(SgNode* node, std::vector<SgNode*>& siblings);

  // Erase the sub-tree with 'node' as its root from the worker list
  // '_scopNodes'.
  void eraseTreeFromScopNodes(SgNode* node);


  std::vector<SgNode*>& getScopRoots();
  std::set<SgNode*> getScopNodes();
  void reinitializeAttributes();
  // Remove all scops which does not contain at least one for loop.
  void pruneScopWithoutLoop();

  // Build/Finalize all necessary control information.
  void
  buildControlInformation();

  /**
   * Accessor.
   *
   */
  std::map<SgNode*, std::vector<SgNode*> > getControlMap();

private:
  std::set<SgNode*> _scopNodes;
  std::vector<SgNode*> _scopRoots;
  std::set<SgNode*> _nonScopNodes;
  SgNode* _rootNode;
  ScopParser _scopParser;
};



#endif

