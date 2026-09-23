from sklearn.metrics import confusion_matrix


def tss(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    tss_sampling = tp / (tp + fn) + tn / (tn + fp) - 1
    return tss_sampling
