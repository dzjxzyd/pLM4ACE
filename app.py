import os
import numpy as np
import pandas
import pickle
from flask import Flask, request, url_for, redirect, render_template, send_from_directory
import joblib
import pandas as pd
from werkzeug.utils import secure_filename

app = Flask(__name__)


# embeddings function
def esm_embeddings(peptide_sequence_list: list):
    # NOTICE: ESM for embeddings is quite RAM usage, if your sequence is too long,
    #         or you have too many sequences for transformation in a single converting,
    #         you conputer might automatically kill the job.
    # return a panda.dataframe
    import torch
    import pandas as pd
    import esm
    import collections
    # load the model
    # NOTICE: if the model was not downloaded in your local environment, it will automatically download it.
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    # load the peptide sequence list into the bach_converter
    batch_labels, batch_strs, batch_tokens = batch_converter(peptide_sequence_list)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    ## batch tokens are the embedding results of the whole data set

    # Extract per-residue representations (on CPU)
    with torch.no_grad():
        # Here we export the last layer of the EMS model output as the representation of the peptides
        # model'esm2_t6_8M_UR50D' only has 6 layers, and therefore repr_layers parameters is equal to 6
        results = model(batch_tokens, repr_layers=[6], return_contacts=True)
    token_representations = results["representations"][6]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
    sequence_representations = []
    for i, tokens_len in enumerate(batch_lens):
        sequence_representations.append(token_representations[i, 1: tokens_len - 1].mean(0))
    # save dataset
    # sequence_representations is a list and each element is a tensor
    embeddings_results = collections.defaultdict(list)
    for i in range(len(sequence_representations)):
        # tensor can be transformed as numpy sequence_representations[0].numpy() or  sequence_representations[0].to_list
        each_seq_rep = sequence_representations[i].tolist()
        for each_element in each_seq_rep:
            embeddings_results[i].append(each_element)
    embeddings_results = pd.DataFrame(embeddings_results).T
    return embeddings_results


# collect the output
def assign_activity(predicted_class):
    import collections
    out_put = []
    for i in range(len(predicted_class)):
        if predicted_class[i] == 0:
            # out_put[int_features[i]].append(1)
            out_put.append('active')
        else:
            # out_put[int_features[i]].append(2)
            out_put.append('non-active')
    return out_put


def get_filetype(filename):
    return filename.rsplit('.', 1)[1].lower()


def model_selection(num: str):
    model = ''
    if num == '1':
        model = 'LR.pkl'
    elif num == '2':
        model = 'SVM.pkl'
    elif num == '3':
        model = 'MLP.pkl'
    return model


def text_fasta_reading(file_name):
    """
    A function for reading txt and fasta files
    """
    import collections
    # read txt file with sequence inside
    file_read = open(file_name, mode='r')
    file_content = []  # create a list for the fasta content temporaty storage
    for line in file_read:
        file_content.append(line.strip())  # extract all the information in the file and delete the /n in the file

    # build a list to collect all the sequence information
    sequence_name_collect = collections.defaultdict(list)
    for i in range(len(file_content)):
        if '>' in file_content[i]:  # check the symbol of the
            sequence_name_collect[file_content[i]].append(file_content[i + 1])

    # transformed into the same style as the xlsx file loaded with pd.read_excel and sequence_list = dataset['sequence']
    sequence_name_collect = pd.DataFrame(sequence_name_collect).T
    sequence_list = sequence_name_collect[0]
    return sequence_list


# create an app object using the Flask class
@app.route('/')
def home():
    return render_template('index.html')


@app.route('/predict', methods=['POST'])
def predict():
    # 每一个网页上的 输入的框，是一个单独的x，下面这个就是吧这个单独的信息变成一个list，每一个单独的就是一个str （也可以吧x变成int 如果想要的话）
    # int_features  = [str(x) for x in request.form.values()] # this command basically use extract all the input into a list
    # final_features = [np.array(int_features)]
    int_features = [str(x) for x in request.form.values()]
    # we have two input in the website, one is the model type and other is the peptide sequences

    # choose scaler and model
    #    name = int_features[0]
    if int(int_features[0]) < 1 or int(int_features[0]) > 12:
        return render_template('index.html')
    model_name = model_selection(int_features[0])
    model=pickle.load(open(model_name,'rb'))
    scaler_name = model_name + '.joblib'
    scaler = joblib.load(scaler_name)

    sequence_list = int_features[1].split(',')  # 因为这个list里又两个element我们需要第二个，所以我只需要把吧这个拿出来，然后split
    # 另外需要注意，这个地方，网页上输入的时候必须要是AAA,CCC,SAS, 这个格式，不同的sequence的区分只能使用逗号，其他的都不可以
    peptide_sequence_list = []
    for seq in sequence_list:
        format_seq = [seq, seq]  # the setting is just following the input format setting in ESM model, [name,sequence]
        tuple_sequence = tuple(format_seq)
        peptide_sequence_list.append(
            tuple_sequence)  # build a summarize list variable including all the sequence information

    embeddings_results = esm_embeddings(peptide_sequence_list)  # conduct the embedding
    normalized_embeddings_results = scaler.transform(embeddings_results)  # normalized the embeddings

    # prediction
    predicted_protability = model.predict(normalized_embeddings_results, batch_size=1)
    predicted_class = []
    for i in range(predicted_protability.shape[0]):
        index = np.where(predicted_protability[i] == np.amax(predicted_protability[i]))[0][0]
        predicted_class.append(index)  # get the class of the results
    predicted_class = assign_activity(predicted_class)  # transform results (0 and 1) into 'active' and 'non-active'
    final_output = []
    for i in range(len(sequence_list)):
        temp_output=sequence_list[i]+': '+predicted_class[i]+';'
        final_output.append(temp_output)

    return render_template('index.html',
                           prediction_text="Prediction results of input sequences {}".format(final_output))


@app.route('/pred_with_file', methods=['POST'])
def pred_with_file():
    # delete existing files that are in the 'input' folder
    dir = 'input'
    for f in os.listdir(os.path.join(os.getcwd(), dir)):
        os.remove(os.path.join(dir, f))
    # 每一个网页上的 输入的框，是一个单独的x，下面这个就是吧这个单独的信息变成一个list，每一个单独的就是一个str （也可以吧x变成int 如果想要的话）
    # int_features  = [str(x) for x in request.form.values()] # this command basically use extract all the input into a list
    # final_features = [np.array(int_features)]
    features = request.form  # .values()
    # we have two input in the website, one is the model type and other is the peptide sequences
    # choose scaler and model
    #    name = int_features[0]
    model_name = model_selection(features.get("Model_selection"))
    model=pickle.load(open(model_name,'rb'))
    scaler_name = model_name + '.joblib'
    scaler = joblib.load(scaler_name)

    file = request.files["Peptide_sequences"]
    filename = secure_filename(file.filename)
    filetype = get_filetype(filename)
    save_location = os.path.join('input', filename)
    file.save(save_location)

    sequence_list = []
    if filetype == 'xls' or filetype == 'xlsx':
        df = pandas.read_excel(save_location, header=0)
        sequence_list = df["sequence"].tolist()
    if filetype == 'txt' or filetype == 'fasta':
        sequence_list = text_fasta_reading(save_location)

    if len(sequence_list) == 0:
        return render_template("index.html")

    # 因为这个list里又两个element我们需要第二个，所以我只需要把吧这个拿出来，然后split
    # 另外需要注意，这个地方，网页上输入的时候必须要是AAA,CCC,SAS, 这个格式，不同的sequence的区分只能使用逗号，其他的都不可以
    peptide_sequence_list = []
    for seq in sequence_list:
        format_seq = [seq, seq]  # the setting is just following the input format setting in ESM model, [name,sequence]
        tuple_sequence = tuple(format_seq)
        peptide_sequence_list.append(
            tuple_sequence)  # build a summarize list variable including all the sequence information

    embeddings_results = esm_embeddings(peptide_sequence_list)  # conduct the embedding
    normalized_embeddings_results = scaler.transform(embeddings_results)  # normalized the embeddings

    # prediction
    predicted_protability = model.predict(normalized_embeddings_results, batch_size=1)
    predicted_class = []
    for i in range(predicted_protability.shape[0]):
        index = np.where(predicted_protability[i] == np.amax(predicted_protability[i]))[0][0]
        predicted_class.append(index)  # get the class of the results
    predicted_class = assign_activity(predicted_class)  # transform results (0 and 1) into 'active' and 'non-active'

    report = {"sequence": sequence_list, "activity": predicted_class}
    report_df = pandas.DataFrame(report)
    save_result_path = os.path.join('input', "report.xlsx")
    report_df.to_excel(save_result_path)
    send_from_directory("input", "report.xlsx")

    return send_from_directory("input", "report.xlsx")


if __name__ == '__main__':
    app.run()
