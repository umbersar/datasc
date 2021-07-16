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
