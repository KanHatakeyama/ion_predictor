# Generated by Django 3.2.12 on 2022-03-16 02:36

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ion_manager', '0006_auto_20220316_0235'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='chemical',
            name='company',
        ),
        migrations.RemoveField(
            model_name='chemical',
            name='made_by',
        ),
    ]