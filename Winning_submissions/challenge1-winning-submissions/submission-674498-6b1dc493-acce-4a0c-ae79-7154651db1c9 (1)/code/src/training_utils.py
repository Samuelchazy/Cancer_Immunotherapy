import numpy as np
import tensorflow as tf
import gc

from src.model import SimpNet
from src.loss import loss_function
from src.eval import evaluate


def make_ds(x, q, y):
    ds = tf.data.Dataset.from_tensor_slices({'q': q, 'x': x, 'y': y})
    ds = ds.shuffle(100_000).repeat()
    return ds


def create_datagen(x, q, y, batch_size, balanced=True):
    strata = [tuple(x) for x in q]
    class_ds = []
    weights = []
    for s in set(strata):
        idx = [s == x for x in strata]
        ds = make_ds(x[idx], q[idx], y[idx])
        class_ds.append(ds)
        weights.append(sum(idx))

    if not balanced:
        weights = np.array(weights)
        weights = weights / weights.sum()
    else:
        weights = np.ones(len(set(strata)))
        weights = weights / weights.sum()

    dataset = tf.data.Dataset.sample_from_datasets(class_ds, weights=weights)
    dataset = dataset.batch(batch_size)
    return dataset.as_numpy_iterator()


def train_model(x_train, q_train, y_train,
                y_test, perturbations_test, condition_test, g2v_embeddings,
                h5_name='temp', batch_size=128, epochs=50,
                learning_rate=1e-5, model_dim=300, pretrained_h5_path=""):
    tf.keras.backend.clear_session()
    model = SimpNet(x_train.shape[-1], model_dim, model_dim)
    optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

    if pretrained_h5_path != "":
        # call to 'create' the model
        model(np.zeros((1, 200), 'float32'),
              np.zeros((1, 64), 'float32'))
        model.load_weights(pretrained_h5_path)

    @tf.function
    def train_step(q, z, x, y):
        with tf.GradientTape() as tape:
            pred_y, pred_x = model(q, z)
            loss = loss_function(y, pred_y, x, pred_x)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        return loss

    @tf.function
    def predict(q, z):
        return model(q, z)

    train_gen = create_datagen(x_train, q_train, y_train, batch_size, balanced=True)

    # training loop
    steps_per_epoch = len(x_train) // batch_size
    best_val = np.inf
    for epoch in range(0, epochs):
        for _ in range(steps_per_epoch):
            batch = next(train_gen)
            batch['z'] = np.random.normal(scale=.5, size=(batch_size, 64))
            batch['z'] = batch['z'].astype('float32')
            train_step(**batch)
        val = evaluate(predict, perturbations_test, condition_test,
                       y_test, g2v_embeddings)
        gc.collect()
        if val < best_val:
            model.save_weights(str(h5_name) + '.h5')
            best_val = val

    return best_val

