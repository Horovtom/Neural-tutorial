// Neural-net-tutorial.cpp : Defines the entry point for the console application.

#include <vector>
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>


class TrainingData {
public:
    TrainingData(const std::string filename);

    bool isEof(void) { return m_trainingDataFile.eof(); }

    void getTopology(std::vector<unsigned> &topology);

    unsigned getNextInputs(std::vector<double> &inputVals);

    unsigned getTargetOutputs(std::vector<double> &targetOutputVals);


private:
    std::ifstream m_trainingDataFile;
};

void TrainingData::getTopology(std::vector<unsigned> &topology) {
    std::string line;
    std::string label;

    std::getline(m_trainingDataFile, line);
    std::stringstream ss(line);
    ss >> label;
    if (this->isEof() || label.compare("topology") != 0) {
        abort();
    }

    while (!ss.eof()) {
        unsigned n;
        ss >> n;
        topology.push_back(n);
    }

    return;
}

TrainingData::TrainingData(const std::string filename) {
    m_trainingDataFile.open(filename.c_str());
}

unsigned TrainingData::getNextInputs(std::vector<double> &inputVals) {
    inputVals.clear();

    std::string line;
    getline(m_trainingDataFile, line);
    std::stringstream ss(line);

    std::string label;
    ss >> label;

    if (label.compare("in:") == 0) {
        double oneValue;
        while (ss >> oneValue) {
            inputVals.push_back(oneValue);
        }
    }

    return (unsigned int) inputVals.size();
}

unsigned TrainingData::getTargetOutputs(std::vector<double> &targetOutputVals) {
    targetOutputVals.clear();

    std::string line;
    getline(m_trainingDataFile, line);
    std::stringstream ss(line);

    std::string label;
    ss >> label;
    if (label.compare("out:") == 0) {
        double oneValue;
        while (ss >> oneValue) {
            targetOutputVals.push_back(oneValue);
        }
    }
    return (unsigned int) targetOutputVals.size();
}

struct Connection {
    double weight, deltaWeight;
};

class Neuron;

typedef std::vector<Neuron> Layer;


//****************** class Neuron ******************
class Neuron {
public:
    Neuron(unsigned numOutputs, unsigned myIndex);

    void feedForward(const Layer &prevLayer);

    void calcOutputGradients(double targetVal);

    void calcHiddenGradients(const Layer &nextLayer);

    void updateInputWeights(Layer &prevLayer);

    void setOutputVal(const double val) { m_outputVal = val; }

    double getOutputVal(void) const { return m_outputVal; }

private:
    static double eta;  //[0,0 ... 1,0] overall net training rate
    static double alpha; //[0,0 ... 1,0] multiplier of last weight change (momentum)
    static double randomWeight(void) { return rand() / double(RAND_MAX); } //This returns random value between 0 and 1
    static double activationFunction(double x);

    static double activationFunctionDerivative(double x);

    double sumDOW(const Layer &nextLayer) const;

    double m_outputVal;
    double m_bias;
    std::vector<Connection> m_outputWeights;
    unsigned m_myIndex;
    double m_gradient;


};

double Neuron::eta = 0.15;
double Neuron::alpha = 0.5;

void Neuron::updateInputWeights(Layer &prevLayer) {
    //The weights to be update are in the Connection container in the neurons in the preceding layer

    for (int i = 0; i < prevLayer.size(); ++i) {
        Neuron &neuron = prevLayer[i];
        double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;
        double newDeltaWeight =
                //Individual input, magnified by the gradient and train rate:
                eta
                * neuron.getOutputVal()
                * m_gradient
                //Also add momentum = a fraction of the previous delta weight
                + alpha
                  * oldDeltaWeight;

        neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
        neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
    }
}

double Neuron::sumDOW(const Layer &nextLayer) const {
    double sum = 0.0;

    for (int i = 0; i < nextLayer.size(); ++i) {
        sum += m_outputWeights[i].weight * nextLayer[i].m_gradient;
    }

    return sum;
}

void Neuron::calcHiddenGradients(const Layer &nextLayer) {
    double dow = sumDOW(nextLayer);
    m_gradient = dow * Neuron::activationFunctionDerivative(m_outputVal);
}

void Neuron::calcOutputGradients(double targetVal) {
    double delta = targetVal - m_outputVal;
    m_gradient = delta * Neuron::activationFunctionDerivative(m_outputVal);
}

//This is gonna be a hyperbolic tangent function
double Neuron::activationFunction(double x) {
    // tanh - output range [-1,0...1,0]
    return tanh(x);
}

double Neuron::activationFunctionDerivative(double x) {
    // tanh derivative
    double tan = tanh(x);
    return 1.0 - tan * tan;
}

void Neuron::feedForward(const Layer &prevLayer) {
    double sum = 0.0;

    //Sum outputs from previous layer and add bias to it
    for (unsigned n = 0; n < prevLayer.size(); ++n) {
        sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
    }

    sum += m_bias;

    //Applying a static activationfunction
    m_outputVal = Neuron::activationFunction(sum);
}

Neuron::Neuron(unsigned numOutputs, unsigned myIndex) {
    m_myIndex = myIndex;
    for (unsigned c = 0; c < numOutputs; ++c) {
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight();

    }
    m_bias = 1;
}

//****************** class Net ******************
class Net {
public:
    Net(const std::vector<unsigned> &topology);

    void feedForward(const std::vector<double> &inputVals);

    void backProp(const std::vector<double> &targetVals);

    //getResults(resultVals) is a function that only reads the result, does not modify the Net object at all. Thats why its const
    void getResults(std::vector<double> &resultVals) const;

    double getRecentAverageError(void) const { return m_recentAverageError; }

private:
    std::vector<Layer> m_layers; //m_layers[layerNum][neuronNum]
    double m_error;
    double m_recentAverageError;
    static double m_recentAverageSmoothingFactor;
};

double Net::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over

void Net::getResults(std::vector<double> &resultVals) const {
    resultVals.clear();

    for (int i = 0; i < m_layers.back().size(); ++i) {
        resultVals.push_back(m_layers.back()[i].getOutputVal());
    }
}

void Net::backProp(const std::vector<double> &targetVals) {
    //Calculate overall net error (RMS - root mean square)
    Layer &outputLayer = m_layers.back();
    m_error = 0.0;

    for (unsigned i = 0; i < outputLayer.size(); ++i) {
        double delta = targetVals[i] - outputLayer[i].getOutputVal();
        m_error += delta * delta;
    }
    m_error /= outputLayer.size();
    m_error = sqrt(m_error);

    //Implement a recent average measurement:

    m_recentAverageError =
            (m_recentAverageError * m_recentAverageSmoothingFactor + m_error) / (m_recentAverageSmoothingFactor + 1.0);

    //Calculate output layer gradients
    //1: The last layer (output layer)
    for (int n = 0; n < outputLayer.size(); ++n) {
        outputLayer[n].calcOutputGradients(targetVals[n]);
    }
    //2: Other layers
    for (unsigned layerNum = (unsigned int) (m_layers.size() - 2); layerNum > 0; --layerNum) {
        Layer &currentLayer = m_layers[layerNum];
        Layer &nextLayer = m_layers[layerNum + 1];

        for (int n = 0; n < currentLayer.size(); ++n) {
            currentLayer[n].calcHiddenGradients(nextLayer);
        }
    }

    //For all layers from outputs to first hidden layer update connection weights
    for (unsigned layerNum = (unsigned int) (m_layers.size() - 1); layerNum > 0; --layerNum) {
        Layer &layer = m_layers[layerNum];
        Layer &prevLayer = m_layers[layerNum - 1];

        for (int n = 0; n < layer.size(); ++n) {
            layer[n].updateInputWeights(prevLayer);
        }
    }

}

void Net::feedForward(const std::vector<double> &inputVals) {
    assert(inputVals.size() == m_layers[0].size());

    //feed the first layer (input neurons)
    for (unsigned i = 0; i < inputVals.size(); ++i) {
        m_layers[0][i].setOutputVal(inputVals[i]);
    }

    //Forward propagate
    for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
        Layer &prevLayer = m_layers[layerNum - 1];
        for (unsigned n = 0; n < m_layers[layerNum].size(); ++n) {
            m_layers[layerNum][n].feedForward(prevLayer);
        }
    }
}

//Define the Net constructor:
Net::Net(const std::vector<unsigned> &topology) {
    unsigned numLayers = (unsigned int) topology.size();
    for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
        //append
        m_layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

        //We made layer: add neurons to layer
        for (unsigned neuronNum = 0; neuronNum < topology[layerNum]; ++neuronNum) {

            //get last element in m_layers (Layer) and append a new neuron to it.
            m_layers.back().push_back(Neuron(numOutputs, neuronNum));
        }

    }
    m_recentAverageError = 0.5;

}


void showVectorVals(std::string label, std::vector<double> &v) {
    std::cout << label << " ";
    for (unsigned i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
    }

    std::cout << std::endl;
}

int main() {
    TrainingData trainData("trainingData.txt");

    std::vector<unsigned> topology;
    trainData.getTopology(topology);
    Net myNet(topology);

    std::vector<double> inputVals, targetVals, resultVals;
    int trainingPass = 0;

    while (!trainData.isEof()) {
        ++trainingPass;
        std::cout << std::endl << "Pass " << trainingPass;

        //Get new input data and feed it forward:
        if (trainData.getNextInputs(inputVals) != topology[0]) break;
        showVectorVals(": Inputs:", inputVals);
        myNet.feedForward(inputVals);

        //Collect the net´s actual output results:
        myNet.getResults(resultVals);
        showVectorVals("Outputs:", resultVals);

        //Train the net what the outputs should have been:
        trainData.getTargetOutputs(targetVals);
        showVectorVals("Targets:", targetVals);
        assert(targetVals.size() == topology.back());

        myNet.backProp(targetVals);

        //Report how well the training is working, average over recent samples:
        std::cout << "Net recent average error: " << myNet.getRecentAverageError() << std::endl;
    }

    std::cout << std::endl << "Done" << std::endl;
}

