import time
import tensorflow as tf
from tensorflow.contrib.rnn import LSTMCell

time_length = 128
batch_size = 400
feature_size = 512
hidden_size = 128

## prepare data
x = tf.random_normal([time_length, batch_size, feature_size], mean=0, stddev=1)
q = tf.FIFOQueue(capacity=4, dtypes=tf.float32)
enqueue_op = q.enqueue(x)
num_threads = 1
qr = tf.train.QueueRunner(q, [enqueue_op] * num_threads)
tf.train.add_queue_runner(qr)
inputs = q.dequeue()
inputs.set_shape(x.get_shape())

y = tf.reduce_mean(tf.reduce_sum(inputs, axis=0), axis=1, keep_dims=True)
labels = tf.cast(tf.greater(y, 0), tf.int32)

## build model
sequence_length = tf.Variable([time_length] * batch_size, dtype=tf.int32)
cell_fw = LSTMCell(num_units=hidden_size)
cell_bw = LSTMCell(num_units=hidden_size)
outputs, state = tf.nn.bidirectional_dynamic_rnn(
    cell_fw=cell_fw,
    cell_bw=cell_bw,
    inputs=inputs,
    sequence_length=sequence_length,
    dtype=tf.float32,
    time_major=True)

outputs_fw, outputs_bw = outputs
outputs = tf.concat([outputs_fw, outputs_bw], axis=2)
outputs = tf.reduce_mean(outputs, axis=0)
outputs = tf.contrib.layers.fully_connected(
    inputs=outputs,
    num_outputs=1,
    activation_fn=None)

losses_op = tf.nn.sigmoid_cross_entropy_with_logits(None, tf.cast(labels, tf.float32), outputs)
losses_op = tf.reduce_mean(losses_op)

y_pred = tf.cast(tf.greater(outputs, 0), tf.int32)
accuracy = tf.reduce_mean(tf.cast(tf.equal(y_pred, labels), tf.float32))
train_op = tf.train.AdamOptimizer(0.001).minimize(losses_op, name="train_op")

t1 = time.time()
with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    coord = tf.train.Coordinator()
    threads = tf.train.start_queue_runners(coord=coord)
    for i in range(50):
        _, losses, acc = sess.run([train_op, losses_op, accuracy])
        print('epoch:%d, loss: %f' % (i, losses))

    coord.request_stop()
    coord.join(threads)
    print("Time taken: %f" % (time.time() - t1))