import pandas as pd


def _get_resolutions_with_single_rank(submission_df_with_accepted_info: pd.DataFrame, rank: str) -> pd.DataFrame:
    """
    Gets records for submissions which only match with the given rank and finds records where the accepted name
    is contained in the submission
    :param submission_df_with_accepted_info:
    :param rank:
    :return:
    """
    submissions_with_rank = submission_df_with_accepted_info[submission_df_with_accepted_info['Accepted_Rank'] == rank]
    submissions_without_rank = submission_df_with_accepted_info[
        submission_df_with_accepted_info['Accepted_Rank'] != rank]

    # Submissions where ranks for all matches is the same
    submissions_with_single_rank = submission_df_with_accepted_info[
        submission_df_with_accepted_info['submitted'].isin(submissions_with_rank['submitted'].values) & ~
        submission_df_with_accepted_info['submitted'].isin(
            submissions_without_rank['submitted'].values)]

    # Matches where accepted name is in submitted name
    submissions_with_single_rank_with_acc_name = submissions_with_single_rank[
        ~submissions_with_single_rank['Accepted_Name'].isna()]
    accepted_names_in_submitted_names = submissions_with_single_rank_with_acc_name[
        submissions_with_single_rank_with_acc_name.swifter.progress_bar(False).apply(lambda x: x['Accepted_Name'] in x['submitted'], axis=1)]

    return accepted_names_in_submitted_names
