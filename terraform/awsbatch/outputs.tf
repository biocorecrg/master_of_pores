output "public_dns" {
  value = aws_instance.classroom.*.public_dns
}

output "instance_id" {
  value = aws_instance.classroom.*.id
}

output "bucket_name" {
  value = aws_s3_bucket.class-bucket.*.bucket
}

output "queue" {
  value = aws_batch_job_queue.spot.name
}

output "queue_gpu" {
  value = aws_batch_job_queue.spot-gpu.name
}

output "rand_string" {
  value = random_string.rand.result
}
