**OTOMO: Outlier-Conscious Transfer-Oriented Minority Oversampling**

OTOMO is an outlier-aware and target-oriented oversampling method for Cross-Project Defect Prediction (CPDP). It generates informative synthetic defective instances by jointly considering source minority structure, target-domain relevance, and outlier information.

Unlike conventional oversampling methods that mainly consider the distribution of the source minority class, OTOMO incorporates information from the unlabeled target project when determining which minority instances should contribute more strongly to synthetic sample generation.

OTOMO is designed around three main principles:

1. **Source minority structure**: OTOMO identifies the minority-class instances in the labeled source project and analyzes their relationships with neighboring minority instances.

2. **Target-oriented relevance**: Information from the unlabeled target project is incorporated into the sampling process so that source minority instances that are more relevant to the target distribution receive greater importance.

3. **Outlier-aware filtering**: Identifies and removes potential outliers before oversampling to prevent atypical instances from influencing synthetic sample generation.
