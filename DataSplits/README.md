# Load Data Splits
The data splits used to train the MPNN, are provided as pickle files. Files terminated in ```class.pkl``` were used in the classification task, while files terminated with ```reg.pkl``` contain the splits used in the regression model. Each file stores a dictionary with three keys:

1. **Train**: Sequences used for training.
2. **Val**: Sequences used to assess the performance at each training step.
3. **Test**: Unseen sequences used to evaluate the performance of the trained model.

```python
import pickle
with open("DataSplit.pkl", "rb") as f:
    data = pickle.load(f)
```
