# model/

`general_pupil_unet.mat` (about 31 MB) goes here. It is **not in git** (`*.mat` is ignored; it is trained on unpublished lab video).

Ways to get it:
- copy it from the lab's shared pipeline folder (`opioid_posthoc_portable\matlab\opioidaddiction-matlab\step01_pupil\one_click\model\`), or
- train it yourself:
  ```
  python -X utf8 inventory_existing_pupil.py
  python -X utf8 make_training_list.py
  ```
  then, in MATLAB, `pupil_train_general('<work>\train_pairs.csv', 'model\general_pupil_unet.mat', 'Epochs', 12, 'TestCsv', '<work>\test_pairs.csv')` (about 40 min on a GPU).
