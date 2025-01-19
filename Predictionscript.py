 #Import the following packages 
    import numpy as np
    import tensorflow as tf
    from tensorflow.keras.models import load_model
    import matplotlib.pyplot as plt
    import pandas as pd
    import pyreadr
    import os
    from sklearn.preprocessing import LabelEncoder
    from joblib import load

#Encoder 
encoder_path = load('Encoder_multinet.joblib')

#Load the three trained models for predicting the stage closest to the target dataset 
    model = load_model('modelP24.keras')
    #model = load_model('modelP48.keras')
    #model = load_model('modelAdult.keras')

#Load average expression values of the training set, for the stage closest to the target dataset
    mean = np.load("meanAdult.npy")
    #mean = np.load("meanP70.npy)
    #mean = np.load("meanP50.npy)
    #mean = np.load("meanP40.npy)
    #mean = np.load("meanP30.npy)
    #mean = np.load("meanP15.npy)

# Import target dataset (see Readme.txt for requirements). These lines are required if you are importing the RDS file. Use the function "pyreadr.read_r" to read the file.
dataset_path = 'Example.rds'  # You can change this.
result = pyreadr.read_r(dataset_path)
Target = result[None]  # Extract the dataframe
Target = Adult.T

# You can alternatively center with the mean of your own dataset. This could provide additional domain adaptation but should only be used if your dataset originates from whole optic lobes, i.e. not a subset achieved by FACS etc.
# meanTarget = Target.mean(axis=0)
# Target -= meanTarget

# Function to make prediction
def predict_with_confidence(model, x, no_classes, n_iter=500):
    prediction = model.predict(x)
    softmax = np.amax(prediction, axis=1)
    classes = np.argmax(prediction, axis=1)

    f = K.function([model.layers[0].input, K.learning_phase()],
               [model.layers[-1].output])

    bootstrap = np.zeros((n_iter,) + (x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f((x, 1))[0]

    BSmax = np.argmax(bootstrap, axis=2)
    confidence = np.sum(BSmax == classes, axis=0) / BSmax.shape[0]
    return classes, softmax, confidence


# If you are classifying a large dataset on a machine with limited memory, use the alternative below instead
def predict_with_confidence(model, x, no_classes, n_iter=100):
    prediction = model.predict(x)
    softmax = np.amax(prediction, axis=1)
    classes = np.argmax(prediction, axis=1)

    @tf.function
    def f(x):
        return model(x, training=True)

    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f(x).numpy()
    BSmax1 = np.argmax(bootstrap, axis=2)

    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f(x).numpy()
    BSmax2 = np.argmax(bootstrap, axis=2)

    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f(x).numpy()
    BSmax3 = np.argmax(bootstrap, axis=2)

    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f(x).numpy()
    BSmax4 = np.argmax(bootstrap, axis=2)

    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))
    for i in range(n_iter):
        bootstrap[i, :, :] = f(x).numpy()
    BSmax5 = np.argmax(bootstrap, axis=2)

    BSmax = np.concatenate((BSmax1, BSmax2, BSmax3, BSmax4, BSmax5), axis=0)
    confidence = np.sum(BSmax == classes, axis=0) / BSmax.shape[0]
    return classes, softmax, confidence

# Predict and save
prediction, softmax, confidence = predict_with_confidence(model, Target, 259)
Targetpreds = encoder.inverse_transform(prediction)
combined_array = np.vstack((Targetpreds, confidence)).T
# Save the combined array to a text file in the specified directory
output_path = os.path.join(directory_path, "TragetPreds_and_Confidence.txt")
np.savetxt(output_path, X=combined_array, delimiter=',', fmt='%s,%.6f', header='Cluster,Confidence', comments='')

# This will optionally plot the confidence values for cells, class-by-class
plt.figure(figsize=(20,10))
plt.scatter(Targetpreds, confidence, alpha=0.1)
plt.xlabel("Class")
plt.ylabel("Confidence")
plt.savefig(os.path.join(directory_path, 'TargetConfidenceadult.png'))
plt.show()
