# mu_hat_conv_model
Version with IR matrix issue

Hi Maria,

So this code is very similar to the last code, but has a few key differences:

1) It makes a continous EEGLAB file for each event seperately and save it (CM_step1). 
2) IT creates an IR response matrix for each event seperately but saves in a structure in same directory as before. (CM_step3). 
3) It runs conv modelling using stimuli and responses and |mu_hat| only. (Currently, although commented out are other versions that include outcome etc) (CM_step4). 

The issue with the IR response matric should be occuring in either step 1 (continous EEGLAB file creation) or step 3 (IR matrix)

But I don't see the error. I have looked over and over again now and cannot find what is wrong.
