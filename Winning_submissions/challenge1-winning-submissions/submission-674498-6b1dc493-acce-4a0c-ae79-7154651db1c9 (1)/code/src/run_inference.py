from loguru import logger
import pandas as pd
from pathlib import Path
import typer
import tensorflow as tf
import numpy as np

from src.g2v import load_embeddings
from src.model import SimpNet
from src.eval import sample_dist


def main(
        test_mode: bool = typer.Option(
            False, help="Predict the test genes, otherwise predict the validation genes"
        ),
        submission_save_dir: Path = typer.Option(
            "./data/processed", help="Predict the test genes, otherwise predict the validation genes"
        ),
        model_dir: Path = typer.Option(
            "./data/processed", help="Directory to save the output model weights"
        ),
        n_models: int = typer.Option(
            5, help="Number of models to use in the model-dir"
        ),
        n_samples: int = typer.Option(
            1000, help="Path to the Gene2Vec embeddings"
        )
):
    genes = ['Aqr', 'Bach2', 'Bhlhe40']
    submission_name = "validation_output.csv"
    if test_mode:
        genes = ['Ets1', 'Fosb', 'Mafk', 'Stat3']
        submission_name = "test_output.csv"

    logger.info(f"Loading embeddings")
    g2v_embeddings = load_embeddings()

    logger.info("Creating model")
    tf.keras.backend.clear_session()
    model = SimpNet(256, 256, 256)

    @tf.function
    def predict(q, z):
        return model.call(q, z)

    # call to 'create' the model
    model(np.zeros((n_samples, 200), 'float32'),
          np.zeros((n_samples, 64), 'float32'))

    logger.info("Predicting labels")
    predictions = np.zeros((len(genes), 5), 'float32')
    for i in range(n_models):
        h5_path = model_dir / f"model-{i}.h5"
        model.load_weights(h5_path)
        for j, g in enumerate(genes):
            dist = sample_dist(predict, g2v_embeddings[g.upper()], n_samples)
            dist = dist / len(genes)
            predictions[j] = predictions[j] + dist

    # normalize the predictions
    for i in range(len(genes)):
        predictions[i] = predictions[i] / predictions[i].sum()

    # generate submission
    columns = ['progenitor', 'effector', 'terminal', 'cycling', 'other']
    submission = {"genes": genes}
    for i, c in enumerate(columns):
        submission[c] = predictions[:, i]
    submission = pd.DataFrame(submission)

    submission.to_csv(submission_save_dir / submission_name, index=False)
    logger.success(f"Submission saved to {submission_save_dir / submission_name}")


if __name__ == "__main__":
    typer.run(main)

