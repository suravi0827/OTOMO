# OTOMO: Outlier-Conscious Transfer-Oriented Minority Oversampling

**OTOMO (Outlier-Conscious Transfer-Oriented Minority Oversampling)** is an outlier-aware and target-oriented oversampling method developed for **Cross-Project Defect Prediction (CPDP)**.

OTOMO generates informative synthetic defective instances by jointly considering the **source minority structure**, **target-domain relevance**, and **outlier information**.

Unlike conventional oversampling methods that mainly rely on the distribution of the source minority class, OTOMO incorporates information from the **unlabeled target project** to determine which source minority instances should contribute more strongly to synthetic sample generation.

---

## Key Features

OTOMO is designed around three main principles:

### 1. Source Minority Structure

OTOMO identifies minority-class instances in the labeled source project and analyzes their relationships with neighboring minority instances. This allows synthetic instances to be generated while considering the local structure of the defective class.

### 2. Target-Oriented Relevance

Information from the unlabeled target project is incorporated into the oversampling process. Source minority instances that are more relevant to the target-domain distribution receive greater importance during synthetic sample generation.

### 3. Outlier-Aware Sampling

OTOMO identifies potential outliers in the source minority class and incorporates this information into the sampling process to reduce the influence of atypical instances during synthetic sample generation.

---
## Requirements

The OTOMO implementation requires:

- **MATLAB R2021a or later**
- **Statistics and Machine Learning Toolbox**

The following MATLAB functions are used:

```text
TreeBagger
fitclinear
perfcurve
kmeans
readmatrix
writetable
table
```

### Required Toolbox Mapping

| Function | Requirement |
|---|---|
| `TreeBagger` | Statistics and Machine Learning Toolbox |
| `fitclinear` | Statistics and Machine Learning Toolbox |
| `perfcurve` | Statistics and Machine Learning Toolbox |
| `kmeans` | Statistics and Machine Learning Toolbox |
| `readmatrix` | MATLAB |
| `writetable` | MATLAB |
| `table` | MATLAB |

No additional MATLAB toolboxes are required for the provided OTOMO implementation.

---

## Running OTOMO

The OTOMO experiments can be executed using the `main.m` script.

### Step 1: Select a Dataset

Open `main.m` and change the dataset name in **Line 5**.

For example:

```matlab
dataset = "NASA";
```

Replace `"NASA"` with the dataset you want to evaluate.

### Step 2: Set the Gamma Value

Change the value of `gamma` in **Line 7** of `main.m`.

For example:

```matlab
gamma = 0.5;
```

Set this value according to the desired experimental configuration.

### Step 3: Run the Experiment

Run:

```matlab
main
```

OTOMO will perform the experiment using the selected dataset and gamma value.

---

## Output

After the experiment is completed, the results are automatically saved as a CSV file.

The output filename follows the format:

```text
dataset_OTOMO_Gamma_value_Results.csv
```

For example:

```text
NASA_OTOMO_Gamma_0.5_Results.csv
```


---

## Repository Structure

```text
OTOMO/
│
├── main.m
├── data/
│   └── ...
│
├── OTOMO-related functions
│   └── ...
│
└── README.md
```

`main.m` is the primary script for configuring and running the experiments.

---

## Citation

If you use OTOMO or this implementation in your research, please cite:

> **OTOMO: Outlier-Conscious Transfer-Oriented Minority Oversampling for Cross-Project Defect Prediction**

The complete bibliographic information will be added after publication.

---



## Contact

For questions regarding OTOMO, its implementation, or the experimental setup, please contact the authors of the paper.
