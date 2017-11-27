Note that the tree building is not done in the Tree constructor; can't
call Tree constructor recursively, which is very natural for trees.
Will this change in new "initializer" approach in 1.17?
