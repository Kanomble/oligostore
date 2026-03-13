from django.shortcuts import get_object_or_404

from ..models import AnalysisJob


def create_analysis_job(*, owner, job_type, sequence_file=None, primer=None, primer_pair=None):
    target_object_id = None
    if primer is not None:
        target_object_id = primer.id
    elif primer_pair is not None:
        target_object_id = primer_pair.id

    return AnalysisJob.objects.create(
        owner=owner,
        job_type=job_type,
        sequence_file=sequence_file,
        primer=primer,
        primer_pair=primer_pair,
        target_object_id=target_object_id,
    )


def mark_job_running(job: AnalysisJob, celery_task_id: str):
    job.status = AnalysisJob.STATUS_RUNNING
    job.celery_task_id = celery_task_id
    job.error_message = ""
    job.save(update_fields=["status", "celery_task_id", "error_message", "updated_at"])
    return job


def mark_job_success(job_id: int, payload):
    job = AnalysisJob.objects.get(id=job_id)
    job.status = AnalysisJob.STATUS_SUCCESS
    job.result_payload = payload
    job.error_message = ""
    job.save(update_fields=["status", "result_payload", "error_message", "updated_at"])
    return job


def mark_job_failure(job_id: int, message: str):
    job = AnalysisJob.objects.get(id=job_id)
    job.status = AnalysisJob.STATUS_FAILURE
    job.error_message = message
    job.save(update_fields=["status", "error_message", "updated_at"])
    return job


def get_owned_job_or_404(*, owner, job_id: int):
    return get_object_or_404(AnalysisJob, id=job_id, owner=owner)
