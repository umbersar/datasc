library(keras)

mnist <- dataset_mnist()
str(mnist)
train_images <- mnist$train$x
train_labels <- mnist$train$y

train_images <- array_reshape(train_images, c(60000, 28 * 28))
hist(train_images,8)
max(train_images)
min(train_images)
train_images <- train_images/255

# train_ds <-cbind.data.frame(train_images, train_labels)
# s1 <- as.data.frame(array_reshape(mnist$train$x, c(60000, 28 * 28)))
# identical(train_ds,s1)
# all.equal(train_ds,s1)

test_images <- mnist$test$x
test_labels <- mnist$test$y

test_images <- array_reshape(test_images, c(10000, 28 * 28))
test_images <- test_images / 255

train_labels <- to_categorical(train_labels)
test_labels <- to_categorical(test_labels)

network <- keras_model_sequential() %>%  
            layer_dense(units = 512, activation = "relu", input_shape = c(28 * 28)) %>%  
            layer_dense(units = 10, activation = "softmax")            

network %>% compile(
  optimizer = "rmsprop",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

# MinMax <- function(x) (x-min(x))/(max(x)-min(x))
# a1<-data.frame(apply(train_ds, 2, MinMax))#this wont work because some columns only have 0 range
# a2<-train_ds/255
# identical(a1,a2)
# all.equal(a1,a2)

network %>% fit(train_images, train_labels, epochs = 5, batch_size = 128)
metrics <- network %>% evaluate(test_images,test_labels)

network %>% predict_classes(test_images[1:10,])

library(keras)
mnist <- dataset_mnist()
train_images <- mnist$train$x
train_labels <- mnist$train$y
test_images <- mnist$test$x
test_labels <- mnist$test$y

dim(train_images)
length(dim(train_images))
str(dim(train_images))
train_images[5,,]
plot(as.raster(train_images[5,,], max = 255))

x <- array(round(runif(1000, 0, 9)), dim = c(64, 3, 32, 10))      
y <- array(5, dim = c(32, 10))                                    
z <- sweep(x, c(3, 4), y, pmax)                                   

library(keras)

imdb <- dataset_imdb(num_words = 10000)
str(imdb)
c(c(train_data, train_labels), c(test_data, test_labels)) %<-% imdb
str(train_data[2])
str(train_labels[[1]])
max(lapply(train_data, max))
max(sapply(train_data, max))

word_index <- dataset_imdb_word_index()                                    
reverse_word_index <- names(word_index)                                    
names(reverse_word_index) <- word_index
decoded_review <- sapply(train_data[[1]], function(index) {                
  word <- if (index >= 3) reverse_word_index[[as.character(index - 3)]]
  if (!is.null(word)) word else "?"
})

#this is one hot encoding
vectorize_sequences <- function(sequences, dimension = 10000) {
  results <- matrix(0, nrow = length(sequences), ncol = dimension)      
  for (i in 1:length(sequences))
    results[i, sequences[[i]]] <- 1                                     
  results
}
vectorize_sequences(list(c(1,2),c(2,3),0,0,0))
#this is one hot encoding
vectorize_sequences_my <- function(sequences, dimension = 10) {
  results <- matrix(0, nrow = length(sequences), ncol = dimension)      
  for (i in 1:length(sequences))
    results[i, sequences[[i]]] <- sequences[[i]]
  results
}
vectorize_sequences_my(list(c(1,2),c(2,3),0,0,0))

#todo:Test it with vectorize_sequences_my encoding
x_train <- vectorize_sequences(train_data)
x_test <- vectorize_sequences(test_data)

y_train <- as.numeric(train_labels)
y_test <- as.numeric(test_labels)

library(keras)

model <- keras_model_sequential() %>%
  layer_dense(units = 16, activation = "relu", input_shape = c(10000)) %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

val_indices <- 1:10000

x_val <- x_train[val_indices,]
partial_x_train <- x_train[-val_indices,]
y_val <- y_train[val_indices]
partial_y_train <- y_train[-val_indices]

model %>% compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

history <- model %>% fit(
  partial_x_train,
  partial_y_train,
  epochs = 20,
  batch_size = 512,
  validation_data = list(x_val, y_val)
)
              
str(history)
plot(history)
results <- model %>% evaluate(x_test, y_test)
results
#re train the model using just 4 epochs to prevent overfitting
model <- keras_model_sequential() %>%
  layer_dense(units = 16, activation = "relu", input_shape = c(10000)) %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")


model %>% compile(
  optimizer = "rmsprop",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

history <- model %>% fit(
  partial_x_train,
  partial_y_train,
  epochs = 4,
  batch_size = 512,
  validation_data = list(x_val, y_val)
)
str(history)
plot(history)

results <- model %>% evaluate(x_test, y_test)
results

model %>% predict(x_test[1:10,])

library(keras)

reuters <- dataset_reuters(num_words = 10000)
c(c(train_data, train_labels), c(test_data, test_labels)) %<-% reuters
length(train_data)
length(test_data)

vectorize_sequences <- function(sequences, dimension = 10000) {
  results <- matrix(0, nrow = length(sequences), ncol = dimension)
  for (i in 1:length(sequences))
    results[i, sequences[[i]]] <- 1
  results
}

x_train <- vectorize_sequences(train_data)            
x_test <- vectorize_sequences(test_data)              

to_one_hot <- function(labels, dimension = 46) {
  results <- matrix(0, nrow = length(labels), ncol = dimension)
  for (i in 1:length(labels))
    results[i, labels[[i]] + 1] <- 1
  results
}

one_hot_train_labels <- to_one_hot(train_labels)          
one_hot_test_labels <- to_one_hot(test_labels)            

one_hot_train_labels_1 <- to_categorical(train_labels)
one_hot_test_labels_1 <- to_categorical(test_labels)

model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = "relu", input_shape = c(10000)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 46, activation = "softmax")

model %>% compile(
  optimizer = "rmsprop",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

val_indices <- 1:1000

x_val <- x_train[val_indices,]
partial_x_train <- x_train[-val_indices,]

y_val <- one_hot_train_labels[val_indices,]
partial_y_train = one_hot_train_labels[-val_indices,]

history <- model %>% fit(
  partial_x_train,
  partial_y_train,
  epochs = 20,
  batch_size = 512,
  validation_data = list(x_val, y_val)
)
plot(history)

#20 epochs result in overfitting..try with 9
model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = "relu", input_shape = c(10000)) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 46, activation = "softmax")

model %>% compile(
  optimizer = "rmsprop",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

history <- model %>% fit(
  partial_x_train,
  partial_y_train,
  epochs = 9,
  batch_size = 512,
  validation_data = list(x_val, y_val)
)

results <- model %>% evaluate(x_test, one_hot_test_labels)
results

test_labels_copy <- test_labels
test_labels_copy <- sample(test_labels_copy)
length(which(test_labels == test_labels_copy)) / length(test_labels)

predictions <- model %>% predict(x_test)
dim(predictions)
predictions[1,]
sum(predictions[1,])
max(predictions[1,])
which(predictions[1,]==max(predictions[1,]))
which.max(predictions[1,])

model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = "relu", input_shape = c(10000)) %>%
  layer_dense(units = 4, activation = "relu") %>%
  layer_dense(units = 46, activation = "softmax")

model %>% compile(
  optimizer = "rmsprop",
  loss = "categorical_crossentropy",
  metrics = c("accuracy")
)

model %>% fit(
  partial_x_train,
  partial_y_train,
  epochs = 20,
  batch_size = 128,
  validation_data = list(x_val, y_val)
)

library(keras)

dataset <- dataset_boston_housing()
c(c(train_data, train_targets), c(test_data, test_targets)) %<-% dataset
str(train_data)
str(test_data)
str(train_targets)

mean <- apply(train_data, 2, mean)                                  
std <- apply(train_data, 2, sd)
train_data <- scale(train_data, center = mean, scale = std)         
test_data <- scale(test_data, center = mean, scale = std)

build_model <- function() {                                1
  model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu",
                input_shape = dim(train_data)[[2]]) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1)
  model %>% compile(
    optimizer = "rmsprop",
    loss = "mse",
    metrics = c("mae")
  )
}
