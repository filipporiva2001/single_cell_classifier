# DATA shuffling
# JupyterLAb (ipynb)

# Adult Shuffle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras.utils import to_categorical
import pyreadr
from joblib import dump
import os

# Define the directory path where your data is located
directory_path = '/n/sci/SCI-004375-NYUDATA/Filippo/Multiome'
os.chdir(directory_path)
print("Current working directory:", os.getcwd())

# Load the RDS file using pyreadr
result = pyreadr.read_r('Adult_markersNN.rds')
df = result[None]  # Extract the pandas DataFrame
data = df.T
print("Data loaded from RDS file:", data.shape)

# Load the CSV file
Clusters = pd.read_csv('Adult_IDs.csv')
Clusters = Clusters['cluster']
print("Clusters loaded from CSV file:", Clusters.shape)

# Shuffle the data and clusters
indices = np.arange(data.shape[0])
np.random.shuffle(indices)
dataShuf = data.values[indices, :]  # Convert to NumPy array and shuffle
ClustersShuf = Clusters.values[indices]
print("Data and Clusters shuffled")
print("Shuffled data shape:", dataShuf.shape)
print("Shuffled clusters shape:", ClustersShuf.shape)

# Save the shuffled indices, data, and clusters
np.save("AdultShufflingIndices.npy", indices)
np.save("AdultData_Shuffled.npy", dataShuf)
np.save("AdultIDs_Shuffled.npy", ClustersShuf)
print("Shuffled data and clusters saved as .npy files")

# Verify the shape of df.columns and the dataShuf shape
print("Number of columns in original DataFrame:", len(df.columns))
print("Number of columns in shuffled data array:", dataShuf.shape[1])

# Ensure the columns match the data shape
if len(df.columns) != dataShuf.shape[1]:
    print("Warning: Column length mismatch. Adjusting column names for shuffled data.")
    columns = [f"Feature_{i}" for i in range(dataShuf.shape[1])]
else:
    columns = df.columns

# Also save the shuffled data and clusters as CSV files
pd.DataFrame(dataShuf, columns=columns).to_csv("AdultData_Shuffled.csv", index=False)
pd.DataFrame(ClustersShuf, columns=['cluster']).to_csv("AdultIDs_Shuffled.csv", index=False)
print("Shuffled data and clusters saved as .csv files")

# Reload the shuffled data and clusters
data = np.load("AdultData_Shuffled.npy")
Clusters = np.load("AdultIDs_Shuffled.npy")
print("Shuffled data and clusters reloaded")

# Encode the clusters
encoder = LabelEncoder()
encoder.fit(Clusters)
Clusters = encoder.transform(Clusters)
Clusters = to_categorical(Clusters)
dump(encoder, 'Encoder.joblib')
print("Clusters encoded and encoder saved")




# Normalize the data
meanData = data.mean(axis=0)
data -= meanData
np.save("AdultMeans.npy", meanData)
print("Data normalized and mean saved")
