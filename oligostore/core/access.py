from django.contrib.auth.models import User
from django.db.models import Q, QuerySet

from .models import PCRProduct, Primer, PrimerPair, Project, SequenceFile


def accessible_primers(user: User) -> QuerySet[Primer]:
    if user is None:
        return Primer.objects.none()
    return Primer.objects.filter(users=user).distinct()


def editable_primers(user: User) -> QuerySet[Primer]:
    if user is None:
        return Primer.objects.none()
    return Primer.objects.filter(creator=user)


def accessible_primer_pairs(user: User) -> QuerySet[PrimerPair]:
    if user is None:
        return PrimerPair.objects.none()
    return PrimerPair.objects.filter(users=user).distinct()


def editable_primer_pairs(user: User) -> QuerySet[PrimerPair]:
    if user is None:
        return PrimerPair.objects.none()
    return PrimerPair.objects.filter(creator=user)


def accessible_projects(user: User) -> QuerySet[Project]:
    if user is None:
        return Project.objects.none()
    return Project.objects.filter(users=user).distinct()


def editable_projects(user: User) -> QuerySet[Project]:
    if user is None:
        return Project.objects.none()
    return Project.objects.filter(creator=user)


def accessible_sequence_files(user: User) -> QuerySet[SequenceFile]:
    if user is None:
        return SequenceFile.objects.none()
    return SequenceFile.objects.filter(Q(uploaded_by=user) | Q(users=user)).distinct()


def editable_sequence_files(user: User) -> QuerySet[SequenceFile]:
    return accessible_sequence_files(user)


def accessible_pcr_products(user: User) -> QuerySet[PCRProduct]:
    if user is None:
        return PCRProduct.objects.none()
    return PCRProduct.objects.filter(Q(creator=user) | Q(users=user)).distinct()


def editable_pcr_products(user: User) -> QuerySet[PCRProduct]:
    if user is None:
        return PCRProduct.objects.none()
    return PCRProduct.objects.filter(creator=user)
