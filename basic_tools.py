import logging
import pandas as pd


def load_data(event_path, info_path):
    logging.info('Data loading ...')
    # Get first events
    df_events = pd.read_csv(event_path)
    # Get demographic information
    df_info = pd.read_csv(info_path)
    df_info['ch_year'] = df_info['date_of_birth'].str.split('-').str[0].astype(float)
    logging.info('Data is loaded.')
    return df_events, df_info

