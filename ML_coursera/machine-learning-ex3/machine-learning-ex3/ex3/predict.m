function p = predict(Theta1, Theta2, X)
%PREDICT Predict the label of an input given a trained neural network
%   p = PREDICT(Theta1, Theta2, X) outputs the predicted label of X given the
%   trained weights of a neural network (Theta1, Theta2)

% Useful values
m = size(X, 1);
num_labels = size(Theta2, 1);

% You need to return the following variables correctly 
p = zeros(size(X, 1), 1);

% ====================== YOUR CODE HERE ======================
% Instructions: Complete the following code to make predictions using
%               your learned neural network. You should set p to a 
%               vector containing labels between 1 to num_labels.
%
% Hint: The max function might come in useful. In particular, the max
%       function can also return the index of the max element, for more
%       information see 'help max'. If your examples are in rows, then, you
%       can use max(A, [], 2) to obtain the max for each row.
%

X = [ones(m, 1) X];
z2 = X * Theta1'
# we are applying the theta for all the different oneVsall classifiers to 
#each example and the resulting matrix is one which have columns=no. of labels and rows=no of 
#examples for which we have to predict

m2 = size(z2, 1);#number of rows in z2 H(z2). H(z2) and z2 are same dimensions
a2 = [ones(m2, 1) sigmoid(z2)];#add dummy feature for calculating next layer nodes/features
  
z3 = a2 * Theta2'
#this will have the prediction probabilities in different cols for a row from a2. taking a max 
#row wise will give us the required p 
a3 = sigmoid(z3)

[max_values,indices] = max(a3, [], 2)
p=indices

% =========================================================================


end
