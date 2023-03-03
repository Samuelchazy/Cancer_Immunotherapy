import tensorflow as tf


def proportion_loss(y, p):
    def get_dist(x):
        x = tf.cast(x, tf.float32)
        x = tf.reduce_sum(x, axis=1)
        x = x / tf.reduce_sum(x)
        return x

    y, p = get_dist(y), get_dist(p)
    loss_ = tf.reduce_sum(tf.abs(y - p))
    return loss_


def loss_function(tar_y, pred_y, tar_x, pred_x):
    loss_ = tf.keras.losses.KLDivergence()(tar_y, pred_y)
    loss_ += tf.keras.losses.mean_absolute_error(tar_x, pred_x)
    return loss_

