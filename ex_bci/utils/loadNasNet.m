function [w1, b1, w2, b2] = loadNasNet(nasNetwork)
% loads the matrixes of a NAS net so that it can be run directly on data
global params
nasNetFilepath = params.nasNetFolderBciComputer;

w1 = load([fullfile(nasNetFilepath, nasNetwork) '_w_hidden']);
b1 = load([fullfile(nasNetFilepath, nasNetwork) '_b_hidden']);
w2 = load([fullfile(nasNetFilepath, nasNetwork) '_w_output']);
b2 = load([fullfile(nasNetFilepath, nasNetwork) '_b_output']);