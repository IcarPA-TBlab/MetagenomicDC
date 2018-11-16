import numpy as np
import sys
np.random.seed(1337)  # for reproducibility
from sklearn.datasets import load_digits
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics.classification import accuracy_score

from dbn import SupervisedDBNClassification

#parameters: sys.argv[1] = input dataset as matrix of k-mers

# Loading dataset
nome_train=sys.argv[1].split(".")[0]				

def load_data(file):
	lista=[]
	records= list(open(file, "r"))
	for seq in records:
		elements=seq.split(",")
		level=elements[-1].split("\n")
		classe=level[0]
		lista.append(classe)

	lista=set(lista)
	classes=list(lista)
	X=[]
	Y=[]
	for seq in records:
		elements=seq.split(",")
		X.append(elements[1:-1])
		level=elements[-1].split("\n")
		classe=level[0]
		Y.append(classes.index(classe))
	X=np.array(X,dtype=float)
	Y=np.array(Y,dtype=int)
	data_max= np.amax(X)
	X = X/data_max
	return X,Y,len(classes),len(X[0])


# Training
def create_model():
    classifier = SupervisedDBNClassification(hidden_layers_structure=[256, 256],
                                             learning_rate_rbm=0.05,
                                             learning_rate=0.1,
                                             n_epochs_rbm=10,
                                             n_iter_backprop=100,
                                             batch_size=32,
                                             activation_function='relu',
                                             dropout_p=0.2,verbose=False)
    return classifier

def train_and_evaluate_model (model, X_train, Y_train, X_test, Y_test):

    model.fit(X_train, Y_train)
    Y_pred = model.predict(X_test)

    
    acc = accuracy_score(Y_test, Y_pred)
    return Y_pred, Y_test


if __name__ == "__main__":
        n_folds = 10
        X_train,Y_train,nb_classes,input_length = load_data(sys.argv[1])
        i=1
        kfold = StratifiedKFold(n_splits=n_folds, shuffle=True)
        for train, test in kfold.split(X_train, Y_train):
                model = None # Clearing the NN.
                model = create_model()
                pred, Y_test = train_and_evaluate_model(model, X_train[train], Y_train[train], X_train[test], Y_train[test])
                np.save("./results/preds_"+nome_train+"_"+str(i),pred)
                np.save("./results/test_"+nome_train+"_"+str(i),Y_test)
                i = i+1






