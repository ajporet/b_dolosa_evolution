library(ape)
tree<-ape::read.tree('/Users/Alexandra/Dropbox (MIT)/Lieberman Lab/Personal lab notebooks/Alex/Fall 2020/Breseq - full tree exploration/whole_tree_nexus.tree')
PatristicDistMatrix<-cophenetic(tree)
PatristicDistMatrix2<-100*PatristicDistMatrix
write.table(PatristicDistMatrix2, file = '/Users/Alexandra/Dropbox (MIT)/Lieberman Lab/Personal lab notebooks/Alex/Fall 2020/Breseq - full tree exploration/whole_tree_pairwise.txt',sep="\t")
