{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "eef1d1a1-23a1-49ce-8edb-85f169113088",
   "metadata": {},
   "outputs": [],
   "source": [
    "#Import the following packages \n",
    "import numpy as np\n",
    "import tensorflow as tf\n",
    "from tensorflow.keras.models import load_model\n",
    "import matplotlib.pyplot as plt\n",
    "import pandas as pd\n",
    "import pyreadr\n",
    "import os\n",
    "from sklearn.preprocessing import LabelEncoder\n",
    "from joblib import load\n",
    "\n",
    "#Encoder \n",
    "encoder_path = load('Encoder_multinet.joblib')\n",
    "\n",
    "#Load the three trained models for predicting the stage closest to the target dataset \n",
    "model = load_model('modelP24.keras')\n",
    "#model = load_model('modelP48.keras')\n",
    "#model = load_model('modelAdult.keras')\n",
    "\n",
    "#Load average expression values of the training set, for the stage closest to the target dataset\n",
    "mean = np.load(\"meanAdult.npy\")\n",
    "#mean = np.load(\"meanP70.npy)\n",
    "#mean = np.load(\"meanP50.npy)\n",
    "#mean = np.load(\"meanP40.npy)\n",
    "#mean = np.load(\"meanP30.npy)\n",
    "#mean = np.load(\"meanP15.npy\n",
    "\n",
    "\n",
    "#Import target dataset (see Readme.txt for requiremets), those line are required if you are importing the rds file, use the finction \"pyreadr.read_r\" to read file \n",
    "dataset_path = 'Example.rds' #you can Change \n",
    "result = pyreadr.read_r(dataset_path)\n",
    "Target = result[None]  # Extract the dataframe\n",
    "Target = Adult.T\n",
    "\n",
    "#You can alternatively center with the mean of your own dataset. This could provide additional domain adaptation but should only be used if your\n",
    "#dataser originates from whole optic lobes, i.e. not a subset achieved by FACS etc.\n",
    "#meanTarget = Target.mean(axis=0)\n",
    "#Target -= meanTarget\n",
    "\n",
    "#Function to make prediction \n",
    "def predict_with_confidence(model, x, no_classes, n_iter=500):\n",
    "    prediction = model.predict(x)\n",
    "    softmax = np.amax(prediction, axis=1)\n",
    "    classes = np.argmax(prediction, axis=1)\n",
    "   \n",
    "    f = K.function([model.layers[0].input, K.learning_phase()],\n",
    "               [model.layers[-1].output])\n",
    "    \n",
    "    bootstrap = np.zeros((n_iter,) + (x.shape[0], no_classes) )\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i,:, :] = f((x, 1))[0]\n",
    "    \n",
    "    BSmax = np.argmax(bootstrap, axis=2)\n",
    "    confidence = np.sum(BSmax==classes, axis=0)/BSmax.shape[0]\n",
    "    return classes, softmax, confidence\n",
    "\n",
    "\n",
    "# If you are classifiying a large dataset on a machine with limited memor, use the alternative below instead \n",
    "def predict_with_confidence(model, x, no_classes, n_iter=100):\n",
    "    prediction = model.predict(x)\n",
    "    softmax = np.amax(prediction, axis=1)\n",
    "    classes = np.argmax(prediction, axis=1)\n",
    "    \n",
    "    @tf.function\n",
    "    def f(x):\n",
    "        return model(x, training=True)\n",
    "    \n",
    "    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i, :, :] = f(x).numpy()\n",
    "    BSmax1 = np.argmax(bootstrap, axis=2)    \n",
    "    \n",
    "    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i, :, :] = f(x).numpy()\n",
    "    BSmax2 = np.argmax(bootstrap, axis=2) \n",
    "    \n",
    "    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i, :, :] = f(x).numpy()\n",
    "    BSmax3 = np.argmax(bootstrap, axis=2) \n",
    "\n",
    "    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i, :, :] = f(x).numpy()\n",
    "    BSmax4 = np.argmax(bootstrap, axis=2) \n",
    "    \n",
    "    bootstrap = np.zeros((n_iter, x.shape[0], no_classes))\n",
    "    for i in range(n_iter):\n",
    "        bootstrap[i, :, :] = f(x).numpy()\n",
    "    BSmax5 = np.argmax(bootstrap, axis=2) \n",
    "    \n",
    "    BSmax = np.concatenate((BSmax1, BSmax2, BSmax3, BSmax4, BSmax5), axis=0)\n",
    "    confidence = np.sum(BSmax == classes, axis=0) / BSmax.shape[0]\n",
    "    return classes, softmax, confidence\n",
    "\n",
    "#Predict and save \n",
    "prediction, softmax, confidence = predict_with_confidence(model, Target, 259)\n",
    "Targetpreds = encoder.inverse_transform(prediction)\n",
    "combined_array = np.vstack((Targetpreds, confidence)).T\n",
    "# Save the combined array to a text file in the specified directory\n",
    "output_path = os.path.join(directory_path, \"TragetPreds_and_Confidence.txt\")\n",
    "np.savetxt(output_path, X=combined_array, delimiter=',', fmt='%s,%.6f', header='Cluster,Confidence', comments='')\n",
    "\n",
    "#This will optionally plot the confidence values for cells, class-by-class\n",
    "plt.figure(figsize=(20,10))\n",
    "plt.scatter(Targetpreds, confidence, alpha=0.1)\n",
    "plt.xlabel(\"Class\")\n",
    "plt.ylabel(\"Confidence\")\n",
    "plt.savefig(os.path.join(directory_path, 'TargetConfidenceadult.png'))\n",
    "plt.show()\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.18"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
