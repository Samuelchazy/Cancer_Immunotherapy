import tensorflow as tf


def encoder_three_layer(d_model, dff):
    return tf.keras.Sequential([
        tf.keras.layers.Dense(2 * dff),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.LayerNormalization(),
        tf.keras.layers.ReLU(),

        tf.keras.layers.Dense(dff),
        tf.keras.layers.BatchNormalization(),
        tf.keras.layers.LayerNormalization(),
        tf.keras.layers.ReLU(),

        tf.keras.layers.Dense(d_model)
    ])


def encoder_two_layer(d_model, dff):
    return tf.keras.Sequential([
        tf.keras.layers.Dense(dff),
        tf.keras.layers.ReLU(),

        tf.keras.layers.Dense(d_model)
    ])


class SimpNet(tf.keras.Model):
    def __init__(self, x_dim, enc_dff=128, cls_dff=128):
        super(SimpNet, self).__init__()

        self.enc = encoder_three_layer(enc_dff, enc_dff)
        self.dec = encoder_two_layer(x_dim, 2 * x_dim)
        self.cls = encoder_two_layer(5, cls_dff)

    def call(self, q, z):
        q = 10. * tf.concat([q, z], -1)
        h = self.enc(q)
        x_out = self.dec(h)
        logits = self.cls(h)
        p_out = tf.keras.activations.softmax(logits, -1)

        return p_out, x_out

