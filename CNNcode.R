CNN.fit <- function(Raw.X.train, Raw.Y.train, Raw.X.test, filters.num = 5, kernels.num = 10){
  
  LOOPBACK = ncol(Raw.X.train)
  NumTrain = length(Raw.Y.train); NumTest = nrow(Raw.X.test)
  
  X.input = array(0.0, dim=c(NumTrain, LOOPBACK))
  Y.input = array(0, dim=NumTrain)
  for(fl in 1:NumTrain){
    
    X.input[fl, ] = Raw.X.train[fl,]
    if(Raw.Y.train[fl]==1)
    Y.input[fl] = Raw.Y.train[fl]
  }
  
  X.test.input = array(0.0, dim=c(NumTest, LOOPBACK))
  for(fl in 1:NumTest){
    
    X.test.input[fl, ] = Raw.X.test[fl,]
    
  }
  
  x_train = array_reshape(X.input, c(dim(X.input), 1))
  x_test = array_reshape(X.test.input, c(dim(X.test.input), 1))
  y_train = to_categorical(Y.input, 2)
  
  model2 <- keras_model_sequential()
  model2 %>% 
    layer_conv_1d(filters=filters.num, kernel_size=kernels.num,  activation = "relu",  input_shape=c(LOOPBACK, 1)) %>%
    #layer_global_max_pooling_1d() %>%
    layer_max_pooling_1d(pool_size = 4) %>%
    layer_flatten() %>% 
    layer_dense(units = kernels.num, activation = 'relu') %>%
    layer_dense(units = 2, activation = 'softmax')
  summary(model2)
  model2 %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  history2 <- model2 %>% fit(
    x_train, y_train, 
    epochs = 60, batch_size = 30, 
    validation_split = 0.2
  )
  predClass <- apply(model2 %>% predict(x_test),1,which.max)
  
  return(list(yhat=(predClass-1)))

}
  