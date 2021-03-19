# Vision

- Juul et al presented a useful <blah blah> and cautioned that standard pointwise averages may not capture ...
- They present a few specific ideas for estimating the 'central set' out of an ensemble of curves
- The idea of the functional boxplot goes back farther ... Juul's method is a special case of *functional depth* (for which fast computational approximations already exist ...). Functional depth is useful because fast and robust. (Why do Juul et al depart from the 'standard' J=2 approach ... ??)
- We can generalize this (possibly at the cost of computational efficiency) to any functional (curvewise) distance, and measure centrality either as the distance to a centroid/Fréchet mean or as average distance to all other members of the ensemble (identical for Euclidean distances)
   - Classic curve distances like Fréchet and DTW ignore phase differences
   - We may want to consider phase vs amplitude differences
   - Mahalanobis distance among features/probes (multivariate generalization of Juul's peak-size example)
- Conclusion: Juul et al make an important point, this can be put in a broader context/we point out some additional choices of centrality measures