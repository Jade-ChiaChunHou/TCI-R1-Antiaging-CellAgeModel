################
# R1 Cell Aging
###############

setwd('/Users/hou/Desktop/tci/R1_antiaging/RNA_antiaging_report')


##################################################
#STEP1: Loading and preprocessing the training set 
##################################################

dataset = read.csv('dataset.csv')
dataset = dataset[, 3:23]
dataset_1 = dataset[, 1:6]
dataset_2 = dataset[, 8:21]

dataset = cbind(dataset_1, dataset_2)




###########################################
#STEP2: Loading and Precrossing the testset
###########################################

testset = read.csv('./testset/191004_R1_cellage.csv')

testset = testset[, 3:23]
testset1 = testset[, 1:6]
testset2 = testset[, 8:21]

testset_w0 = cbind(testset1, testset2)



row = seq(1, dim(testset_w0)[1], 3)
avg_testset_w0 = matrix(NA, nrow = length(row), ncol = 20)


for(j in 1:length(row)){
  for (i in 1:20){
    
    avg_testset_w0[j,i] = mean(testset_w0[row[j]:(row[j]), i], na.rm = T)
    
  }
}

avg_testset_w0[,20] = 22

nor_testset_w0 = avg_testset_w0[,1:19] / avg_testset_w0[,20]
nor_testset_w0 = data.frame(nor_testset_w0)
nor_testset_w0 = t(nor_testset_w0)
colnames(nor_testset_w0) = colnames(testset_w0)[1:19]

colnames(dataset) == colnames(nor_testset_w0)
colnames(dataset)
colnames(nor_testset_w0)

testset_1 = nor_testset_w0[1:9]
testset_2 = nor_testset_w0[11]
testset_3 = nor_testset_w0[10]
testset_4 = nor_testset_w0[12:16]
testset_5 = nor_testset_w0[19]
testset_6 = nor_testset_w0[17:18]

final_testset_w0 = c(testset_1, testset_2, testset_3, testset_4, testset_5, testset_6)

colnames(dataset)
colnames(final_testset_w0)
final_testset_w0 = data.frame(final_testset_w0)
final_testset_w0 = t(final_testset_w0)
colnames(final_testset_w0) = colnames(dataset)[1:19]



##################################################
#STEP3: Predict the cell age (Model: RandomForest)
##################################################

# enter the true age
true_age = 73


library(randomForest)
set.seed(1234)
regressor = randomForest(x = dataset[1:19],
                         y = dataset$age,
                         ntree = 500)
regressor
summary(regressor)

pred_191004 = predict(regressor, newdata = final_testset_w0)
pred_191004

pred_191004 = pred_191004 + 20
delta =  pred_191004 - true_age


#################################################
#STEP4: gene group 1: antiaging -> predict the cell age
#################################################

library(randomForest)
set.seed(1234)
regressor_g1 = randomForest(x = dataset[1:5],
                            y = dataset$age,
                            ntree = 500)

pred_g1 = predict(regressor_g1, newdata = final_testset_w0[1:5])
pred_g1

pred_g1 = pred_g1 + 20
delta_g1 = pred_g1 - true_age



####################################################
#STEP5: gene group 2: mitochondria -> predict the cell age
####################################################

library(randomForest)
set.seed(1234)
regressor_g2 = randomForest(x = dataset[6:16],
                            y = dataset$age,
                            ntree = 500)

pred_g2 = predict(regressor_g2, newdata = final_testset_w0[6:16])
pred_g2
pred_g2 = pred_g2 + 20

delta_g2 = pred_g2 - true_age



##################################################
#STEP6: gene group 3: telomerase -> predict the cell age
##################################################

library(randomForest)
set.seed(1234)
regressor_g3 = randomForest(x = dataset[17:19],
                            y = dataset$age,
                            ntree = 500)

pred_g3 = predict(regressor_g3, newdata = final_testset_w0[17:19])
pred_g3
pred_g3 = pred_g3 + 20

delta_g3 = pred_g3 - true_age


##############################################################
#STEP7: The Level of 3 gene groups(IDEAL, STANDARD, NON-IDEAL)
##############################################################


# light1: antiaging
light1 = matrix(NA, nrow = 1, ncol = 1)

for(i in 1:1){
  
  if (delta_g1[i] > 3) { 
    light1[i] = "NON-IDEAL"
    
  } else if (delta_g1[i] <= 3 & delta_g1[i] >= -2) {
    light1[i] = "STANDARD"
    
  } else if  (delta_g1[i] < -2) {
    light1[i] = "IDEAL"
  } 
  
}


# light2: mitochondria
light2 = matrix(NA, nrow = 1, ncol = 1)

for(i in 1:1){
  
  if (delta_g2[i] > 3) { 
    light2[i] = "NON-IDEAL"
    
  } else if (delta_g2[i] <= 3 & delta_g2[i] >= -2) {
    light2[i] = "STANDARD"
    
  } else if  (delta_g2[i] < -2) {
    light2[i] = "IDEAL"
  } 
  
}


# light3: telomerase
light3 = matrix(NA, nrow = 1, ncol = 1)

for(i in 1:1){
  
  if (delta_g3[i] > 3) { 
    light3[i] = "NON-IDEAL"
    
  } else if (delta_g3[i] <= 3 & delta_g3[i] >= -2) {
    light3[i] = "STANDARD"
    
  } else if  (delta_g3[i] < -2) {
    light3[i] = "IDEAL"
  } 
  
}