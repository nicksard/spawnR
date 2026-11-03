package com.nicksard.appsleep;

import android.service.notification.NotificationListenerService;
import android.service.notification.StatusBarNotification;

import java.time.LocalDateTime;

/** Optional: dismisses notifications from apps that are asleep (if the schedule asks for it). */
public class NotificationSilencer extends NotificationListenerService {

    @Override
    public void onListenerConnected() {
        try {
            StatusBarNotification[] active = getActiveNotifications();
            if (active != null) for (StatusBarNotification sbn : active) maybeSilence(sbn);
        } catch (SecurityException ignored) {
            // Access was revoked between connecting and asking.
        }
    }

    @Override
    public void onNotificationPosted(StatusBarNotification sbn) {
        maybeSilence(sbn);
    }

    private void maybeSilence(StatusBarNotification sbn) {
        if (sbn == null || !sbn.isClearable()) return;
        Schedule s = ScheduleStore.sleepingScheduleFor(this, sbn.getPackageName(), LocalDateTime.now());
        if (s != null && s.silenceNotifications) cancelNotification(sbn.getKey());
    }
}
