
import os
import pandas as pd
from sklearn.model_selection import LeaveOneOut, train_test_split, StratifiedKFold, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import *
from sklearn.metrics import confusion_matrix
import numpy as np
from numpy import interp
import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

GLOBAL_RANDOM_STATE = 32
n_folds = 5
skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=GLOBAL_RANDOM_STATE)




def loadFeature(Case_path, Control_path, X_train, X_test, Feature):

    def clean_filename(fname, feature_name):

        name = os.path.splitext(fname)[0]

        if feature_name.lower() in name.lower():

            name = name.replace(f"_{feature_name}", "")
            name = name.replace(f"_{feature_name.lower()}", "")

            name = name.replace(feature_name, "")
            name = name.replace(feature_name.lower(), "")


        name = name.strip("_. ")
        return name


    def build_file_map(directory, feature_name):
        file_map = {}
        if not os.path.exists(directory):
            return file_map

        files = os.listdir(directory)
        for f in files:
            if not f.endswith(('.txt', '.csv', '.tsv')): continue

            cleaned_id = clean_filename(f, feature_name)
            file_map[cleaned_id] = f

            raw_base = os.path.splitext(f)[0]
            if raw_base not in file_map:
                file_map[raw_base] = f

        return file_map


    control_map = build_file_map(os.path.join(Control_path, Feature), Feature)
    case_map = build_file_map(os.path.join(Case_path, Feature), Feature)

    train_X, train_y = [], []
    test_X, test_y = [], []


    def read_and_process(file_path):
        try:
            tmp_data = pd.read_csv(file_path, sep='\t', header=None, names=[Feature, 'value'])
            if Feature in ['BPM', 'EDM', 'SBM', 'JNM', 'OJM']:
                f1 = tmp_data['value'].tolist()
                std_val = np.std(f1)
                if std_val == 0:
                    return (f1 - np.mean(f1))
                else:
                    return (f1 - np.mean(f1)) / std_val
            elif Feature in ['CNV_onco', 'CNV_im']:
                return tmp_data['value'].fillna(0).to_list()
            else:
                return tmp_data['value'].tolist()
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
            return None


    for sample_name in X_train:
        found = False
        data = None
        label = -1


        target_file = None


        if sample_name in control_map:
            target_file = control_map[sample_name]
        else:

            for k, v in control_map.items():
                if k.startswith(sample_name) or sample_name in k:
                    target_file = v
                    break

        if target_file:
            path = os.path.join(Control_path, Feature, target_file)
            data = read_and_process(path)
            if data is not None:
                train_X.append(data)
                train_y.append(0)
                found = True


        if not found:
            target_file = None
            if sample_name in case_map:
                target_file = case_map[sample_name]
            else:
                for k, v in case_map.items():
                    if k.startswith(sample_name) or sample_name in k:
                        target_file = v
                        break

            if target_file:
                path = os.path.join(Case_path, Feature, target_file)
                data = read_and_process(path)
                if data is not None:
                    train_X.append(data)
                    train_y.append(1)
                    found = True

        if not found:
            print(f"Warning: Sample '{sample_name}' missing in feature '{Feature}'")

            pass


    for sample_name in X_test:
        found = False
        data = None

        # 找 Control
        target_file = None
        if sample_name in control_map:
            target_file = control_map[sample_name]
        else:
            for k, v in control_map.items():
                if k.startswith(sample_name) or sample_name in k:
                    target_file = v
                    break

        if target_file:
            data = read_and_process(os.path.join(Control_path, Feature, target_file))
            if data is not None:
                test_X.append(data)
                test_y.append(0)
                found = True

        if not found:
            target_file = None
            if sample_name in case_map:
                target_file = case_map[sample_name]
            else:
                for k, v in case_map.items():
                    if k.startswith(sample_name) or sample_name in k:
                        target_file = v
                        break

            if target_file:
                data = read_and_process(os.path.join(Case_path, Feature, target_file))
                if data is not None:
                    test_X.append(data)
                    test_y.append(1)
                    found = True

    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def GridSearch_base(X_train, y_train):
    rf = RandomForestClassifier(random_state=32)
    param_grid = {

        'n_estimators': [10, 20, 30, 50, 100, 200],
        'max_depth': [None, 1, 2, 5, 10, 20, 30],
        'min_samples_split': [2, 5, 7, 10],
        'min_samples_leaf': [1, 2, 4, 6, 8],
        'max_features': ['sqrt', 'log2', None]
    }

    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        cv=skf,
        scoring='accuracy',
        n_jobs=-1,
        verbose=1
    )
    grid_search.fit(X_train, y_train)
    return grid_search.best_params_, grid_search.best_score_, grid_search.best_estimator_


def GridSearch_meta(X_train, y_train):
    meta_param_grid = {
        'n_estimators': [50, 100, 200, 300],
        'max_depth': [5, 10, 15, 20, None],
        'min_samples_split': [2, 5, 7],
        'min_samples_leaf': [1, 2, 3, 4]
    }
    meta_clf = RandomForestClassifier(random_state=GLOBAL_RANDOM_STATE)
    meta_grid_search = GridSearchCV(
        estimator=meta_clf,
        param_grid=meta_param_grid,
        cv=skf,
        scoring='accuracy',
        n_jobs=-1,
        verbose=1
    )
    meta_grid_search.fit(X_train, y_train)
    return meta_grid_search.best_params_, meta_grid_search.best_score_, meta_grid_search.best_estimator_


def training_BaseModel(X_train, y_train, disease_test, best_params, blend_train, j):
    clf = RandomForestClassifier(random_state=32, **best_params)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=GLOBAL_RANDOM_STATE)
    blend_disease_j = np.zeros((disease_test.shape[0], n_folds))

    for i, (train_index, cv_index) in enumerate(skf.split(X_train, y_train)):
        print(f'Fold [{i}] Feature [{j}]')
        tr_X = X_train[train_index]
        tr_y = y_train[train_index]
        cv_X = X_train[cv_index]
        cv_y = y_train[cv_index]
        clf.fit(tr_X, tr_y)
        blend_train[cv_index, j] = clf.predict_proba(cv_X)[:, 1]
        blend_disease_j[:, i] = clf.predict_proba(disease_test)[:, 1]

    return blend_train, blend_disease_j



def main():

    parser = argparse.ArgumentParser(description='Building Classifiers and Testing Classifiers')
    parser.add_argument('Control_samples_dir', help='directory of control sample')
    parser.add_argument('DiseaseCase_samples_dir', help='directory of case sample')
    #parser.add_argument('-o', '--output', default='output.txt',
     #                   help='path to output')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='显示详细输出')

    args = parser.parse_args()

    if args.verbose:
        print(f"Control directory: {args.Control_samples_dir}")
        print(f"Case directory: {args.DiseaseCase_samples_dir}")

    feature_type = ['BPM', 'EDM', 'SBM', 'JNM', 'OJM', 'CNV_onco', 'CNV_im', 'OLR']

    # check path 
    for i in feature_type:
        Control_feature_dir = os.path.join(args.Control_samples_dir, i)
        if not (os.path.exists(Control_feature_dir) and os.path.isdir(Control_feature_dir)):
            print("Error: cannot find the feature directory :" + Control_feature_dir)
            sys.exit(1)
        DiseaseCase_feature_dir = os.path.join(args.DiseaseCase_samples_dir, i)
        if not (os.path.exists(DiseaseCase_feature_dir) and os.path.isdir(DiseaseCase_feature_dir)):
            print("Error: cannot find the feature directory :" + DiseaseCase_feature_dir)
            sys.exit(1)


    base_feature = feature_type[2]
    Control_base_dir = os.path.join(args.Control_samples_dir, base_feature)
    Case_base_dir = os.path.join(args.DiseaseCase_samples_dir, base_feature)

    if not os.path.exists(Control_base_dir) or not os.path.exists(Case_base_dir):
        print(f"Error: Base feature directory for {base_feature} does not exist.")
        sys.exit(1)

    Control_files = os.listdir(Control_base_dir)
    DiseaseCase_files = os.listdir(Case_base_dir)

    all_name = []
    n_label = []


    def clean_name_simple(fname, feat):
        n = os.path.splitext(fname)[0]
        n = n.replace(f"_{feat}", "").replace(feat, "").strip("_. ")
        return n

    for i in Control_files:
        if not i.endswith(('.txt', '.csv', '.tsv')): continue
        sample_name = clean_name_simple(i, base_feature)
        n_label.append(0)
        all_name.append(sample_name)

    for i in DiseaseCase_files:
        if not i.endswith(('.txt', '.csv', '.tsv')): continue
        sample_name = clean_name_simple(i, base_feature)
        n_label.append(1)
        all_name.append(sample_name)

    if len(all_name) == 0:
        print("Error: No samples found. Check file names and feature suffix.")
        sys.exit(1)

    print(f"Total samples identified: {len(all_name)}")

    X_train_names, X_test_names, _, _ = train_test_split(all_name, n_label, stratify=n_label, test_size=0.3,
                                                         random_state=2)


    feature_X_train = {}
    feature_y_train = {}
    feature_X_test = {}
    feature_y_test = {}
    feature_best_params = {}

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)

    num_base_models = len(feature_type)
    blend_train = np.zeros((len(X_train_names), num_base_models))
    blend_test = np.zeros((len(X_test_names), num_base_models))

    j = 0
    roc = []
    results_disease = []
    all_results = []


    for feature in feature_type:
        print(f"\nProcessing feature: {feature}...")


        key_X_tr = f"{feature}_X_train"
        key_y_tr = f"{feature}_y_train"
        key_X_te = f"{feature}_X_test"
        key_y_te = f"{feature}_y_test"


        tmp_X_tr, tmp_y_tr, tmp_X_te, tmp_y_te = loadFeature(
            args.DiseaseCase_samples_dir, args.Control_samples_dir, X_train_names, X_test_names, feature
        )


        if len(tmp_X_tr) == 0:
            print(f"WARNING: Feature '{feature}' has 0 samples. Skipping.")

            j += 1
            continue


        feature_X_train[key_X_tr] = tmp_X_tr
        feature_y_train[key_y_tr] = tmp_y_tr
        feature_X_test[key_X_te] = tmp_X_te
        feature_y_test[key_y_te] = tmp_y_te

        best_params, best_score_cv, best_rf = GridSearch_base(feature_X_train[key_X_tr], feature_y_train[key_y_tr])
        feature_best_params[key_X_tr] = best_params

        blend_train, blend_test_j = training_BaseModel(
            feature_X_train[key_X_tr], feature_y_train[key_y_tr], feature_X_test[key_X_te], best_params, blend_train, j
        )
        blend_test[:, j] = blend_test_j.mean(1)
        j += 1


        y_pred = best_rf.predict(feature_X_test[key_X_te])
        clf_tn, clf_fp, clf_fn, clf_tp = confusion_matrix(feature_y_test[key_y_te], y_pred).ravel()
        accuracy = accuracy_score(feature_y_test[key_y_te], y_pred)
        recall = recall_score(feature_y_test[key_y_te], y_pred)
        specificity = clf_tn / (clf_tn + clf_fp) if (clf_tn + clf_fp) > 0 else 0
        f1 = f1_score(feature_y_test[key_y_te], y_pred)
        clf_prob = best_rf.predict_proba(feature_X_test[key_X_te])[:, 1]

        try:
            clf_roc_auc = roc_auc_score(feature_y_test[key_y_te], clf_prob)
        except ValueError:
            clf_roc_auc = 0.5

        clf_fpr, clf_tpr, _ = roc_curve(feature_y_test[key_y_te], clf_prob)
        clf_mean_fpr1 = np.linspace(0, 1, 100)
        clf_tpr_interp = interp(clf_mean_fpr1, clf_fpr, clf_tpr)
        clf_tpr_interp[0] = 0.0
        clf_tpr_interp[-1] = 1.0

        results_disease.append(
            [feature, args.DiseaseCase_samples_dir, accuracy, specificity, recall, f1, clf_roc_auc]
        )
        roc.append([feature, clf_mean_fpr1, clf_tpr_interp, clf_roc_auc])

        all_results.append({
            'Dataset': args.DiseaseCase_samples_dir,
            'Feature': feature,
            'CV_Accuracy': best_score_cv,
            'Test_Accuracy': accuracy,
            'Best_Parameters': best_params
        })

    print("\nTraining Meta-Model (CFECC)...")


    target_feat = feature_type[2]
    target_key_y_train = f"{target_feat}_y_train"
    target_key_y_test = f"{target_feat}_y_test"


    if target_key_y_train not in feature_y_train:
        print(f"Warning: Target feature '{target_feat}' data missing. Searching for alternative...")
        for f in feature_type:
            if f"{f}_y_train" in feature_y_train:
                target_feat = f
                target_key_y_train = f"{f}_y_train"
                target_key_y_test = f"{f}_y_test"
                break

    print(f"Using labels from feature '{target_feat}' as Ground Truth.")

    y_train_final = feature_y_train[target_key_y_train]
    y_test_final = feature_y_test[target_key_y_test]

    best_meta_params, best_meta_score_cv, best_meta_clf = GridSearch_meta(blend_train, y_train_final)

    blend_pred = best_meta_clf.predict(blend_test)
    blend_prob = best_meta_clf.predict_proba(blend_test)[:, 1]

    clf_tn, clf_fp, clf_fn, clf_tp = confusion_matrix(y_test_final, blend_pred).ravel()
    blend_accuracy = accuracy_score(y_test_final, blend_pred)
    blend_recall = recall_score(y_test_final, blend_pred)
    blend_specificity = clf_tn / (clf_tn + clf_fp) if (clf_tn + clf_fp) > 0 else 0
    blend_f1 = f1_score(y_test_final, blend_pred)

    try:
        blend_auc = roc_auc_score(y_test_final, blend_prob)
    except:
        blend_auc = 0.5

    blend_fpr, blend_tpr, _ = roc_curve(y_test_final, blend_prob)
    blend_mean_fpr1 = np.linspace(0, 1, 100)
    blend_tpr_interp = interp(blend_mean_fpr1, blend_fpr, blend_tpr)
    blend_tpr_interp[0] = 0.0
    blend_tpr_interp[-1] = 1.0

    all_results.append({
        'Dataset': args.DiseaseCase_samples_dir,
        'Feature': 'CFECC',
        'CV_Accuracy': best_meta_score_cv,
        'Test_Accuracy': blend_accuracy,
        'Best_Parameters': best_meta_params
    })
    df_results = pd.DataFrame(all_results)
    df_results.to_csv('GridSearch_All_Features_Results.csv', index=False, encoding='utf-8')

    results_disease.append(
        ['CFECC', args.DiseaseCase_samples_dir, blend_accuracy, blend_specificity, blend_recall, blend_f1, blend_auc]
    )
    roc.append(['CFECC', blend_mean_fpr1, blend_tpr_interp, blend_auc])
    roc.to_csv('Roc.csv', index=False, encoding='utf-8')

    df_disease = pd.DataFrame(results_disease,
                              columns=['Classifier', 'DiseaseCase', 'Accuracy', 'Specificity', 'Recall', 'F1', 'AUC'])
    df_disease.to_csv('ClassifyPerformance.csv', index=False, encoding='utf-8')                 
    print("\n=== Final Results ===")
    print(df_disease)


if __name__ == "__main__":
    main()